"""Extract topological and hydraulic information from EPANET2 models.

Module: epanet
Author: Will Furnass
Date:   2010-11-12



Edited for Node data and plotting 
Richard Collins
15-11-2010

"""
from ctypes import CDLL, byref, c_long, c_int, c_float, create_string_buffer
import sys, os.path

import pylab as pp
import numpy as np

global lib_epanet
class ENNode(object):
	"""A node in an EPANET2 model."""
	def __init__(self, node_name, point):
		"""Create a new node.

		Instance variables:
		node_name -- name/ID of node
		point -- coordinates of node expressed as a list i.e. [x,y] 

		"""
		self.node_name, self.point = node_name, point
		self.head = [0]
		self.monte_head = []
		self.monte_demand = []

	def wkt(self):
		"""Get a (GIS) Well-Known Text representation (as a string)."""
		return "POINT(%f %f)" % (self.point[0],self.point[1])
        
class ENPipe(object):
	"""A pipe in an EPANET2 model."""
	def __init__(self, pipe_name, node_1, node_2, length, diameter, roughness,wall_thickness = 10.,modulus = 2e9):
		"""Create a new pipe.

		Instance variables:
		pipe_name -- name/ID of pipe
		node_1 -- name/ID of one of two terminating nodes
		node_2 -- name/ID of other terminating node
		sim_velocities -- all pipe velocity vectors associated with a 24h extended simulation (as a list of real numbers)
		points -- ordered list of all coodinates (terminating and intermediary); coords expressed in [x,y] form.

		"""   
		self.pipe_name, self.node_1, self.node_2 = pipe_name, node_1, node_2
		self.length, self.diameter, self.roughness, self.wall_thickness, self.modulus = float(length), float(diameter), float(roughness), float(wall_thickness), float(modulus)
		
		self.sim_velocities = []
		self.points = []
		
		####		Calculation of the pipe wave speed (celerity) for transient calculations
		k = 2.05e9
		rho = 998.2
		self.celerity = np.sqrt(1./(rho*(1./k + self.diameter/(self.wall_thickness*self.modulus))))
		self.celerity = 100.		###  Use for VerySimple.inp models

	def add_point(self, p):
		"""Add a new point in [x,y] form to the 'points' instance variable.""" 
		self.points.append(p)

	def mean_velocity(self):
		"""Returns the average diurnal velocity as the mean of the velocity vectors stored in the 'sim_velocities' instance var."""  
		return float(sum(self.sim_velocities))/len(self.sim_velocities)

	def wkt(self):
		"""Get a (GIS) Well-Known Text representation (as a string)."""
		# Returns something like "LINESTRING(0 0,1 1,1 2)"
		return "LINESTRING(%s)" % ','.join("%f %f" % (p[0], p[1]) for p in self.points)

class ENValve(object):
	"""Represents a valve in an EPANET model."""
	def __init__(self, valve_name, node_1, node_2, setting, fixed_setting):
		"""Create a new valve.
		
		Instance variables:
		valve_name -- name/ID of valve
		node_1 -- name/ID of one of two terminating nodes
		node_2 -- name/ID of other terminating node
		setting -- valve setting
		fixed_setting -- valve fixed setting
		"""
		self.valve_name, self.node_1, self.node_2, self.setting, self.fixed_setting = \
		valve_name, node_1, node_2, setting, fixed_setting

class ENModel(object):
    """Represents an EPANET2 model."""
    # The EPANET2 Toolkit DLL is located in the same dir as this module.  
    dll_loc = '/lib/libepanet2.so'  ##Location of the compiled epanet library
    def __init__(self, inp_filename):
        """Create a representation of an EPANET2 model from an EPANET2 input (*.INP) file.
        Davbe
        Arguments:
        inp_filename -- EPANET2 input file (*.INP)
        
        Automatically creates three lists as instance parameters for storing objects that correspond 
        to network assets:    
        nodes  -- list of Node objects which store node names (IDs) and coordinates 
        pipes  -- list of Pipe objects which store pipe names, terminating node pairs, velocities corresponding
                  to each timestep of a 24h extended simulation and pipe geometries
        valves -- list of Valve objects which store valve names, valve terminating node pairs and valve types
        
        The Node and Pipe objects have methods for accessing Well-Known Text GIS representations of their geometries
        and the latter has a method for determining the mean diurnal velocity of a pipe.
        
        """
        
        
        self.inp_filename = inp_filename
        self.nodes, self.pipes, self.valves, self.junctions = [], [], [], []
        # Primary-key-like objects for getting references to ENPipe and ENNode using asset IDs.
        # Used for updating per-pipe lists of vertices  
        self.pipe_idx = {} 
        self.node_idx = {}
        self.junction_idx = {}
        self.pipe_lengths = []
        self._load_static_inp_info()
        #self._load_simulation_results()
        
        
        
    def _load_static_inp_info(self):
        inp_file = open(self.inp_filename, 'r')
        self.lib_epanet = CDLL(ENModel.dll_loc)
        rep_filename = self.inp_filename.lower().rstrip('.inp') + '.rep'
	out_filename = self.inp_filename.lower().rstrip('.inp') + '.out'
	self.lib_epanet.ENopen(self.inp_filename, rep_filename, out_filename)
        self.lib_epanet.ENopenH()
        self.lib_epanet.ENsettimeparam(0, 0)
        self.lib_epanet.ENinitH(0)
        
        time, time_step = c_long(), c_long()
       	i, num_pipes, num_nodes = c_int(), c_int(), c_int()
       	flow, abs_velocity= c_float(), c_float()
       	head,demand = c_float(), c_float()
       	pipe_id = create_string_buffer(32)
       	node_id = create_string_buffer(32)
       	real_pipe_id = ""
       	velocity = 0
        
        self.lib_epanet.ENrunH(byref(time))	
        
        section = None
        for line in inp_file:
            # Ignore rows that are either blank or contain column headings
            if line.startswith(';') or line.strip() == "":
                continue
            # Detect start of section
            elif line.startswith('['):
                section = line.strip()
            else:
                vals = line.split()
                if section == '[PIPES]': # pipe ends (node 1 & 2)
                    new_pipe = ENPipe(vals[0], vals[1], vals[2], vals[3], vals[4], vals[5])
                    self.pipes.append(new_pipe)
                    self.pipe_idx[new_pipe.pipe_name] = new_pipe
                    self.pipe_lengths.append(float(vals[3]))
                    
                elif section == '[VALVES]':
                    self.valves.append(ENValve(vals[0], vals[1], vals[2], vals[5], vals[7]))
                    
                elif section == '[COORDINATES]': # node coords	note this is the junction coordinates, useful for plotting
                    new_node = ENNode(vals[0], [float(vals[1]), float(vals[2])])
                    self.nodes.append(new_node)
                    self.node_idx[new_node.node_name] = new_node 
                    
                elif section == '[VERTICES]': # pipe point coords
                    pipe_name, x, y = vals[0:3]
                    self.pipe_idx[pipe_name].points.append([float(x),float(y)])
                    
#                elif section == '[JUNCTIONS]': #Work out which of the nodes are junctions
#                    new_junction = ENJunction(vals[0],vals[2])
#                    self.junctions.append(new_junction)
#                    self.junction_idx[new_junction.junc_name] = new_junction
                    #print vals
        # Pipe coords in [VERTICES] section don't include coords of terminating nodes 
        # so must prepend coords of node_1 and postpend coords of node_2 to pipe coord lists. 
        for pipe in self.pipes:
            pipe.points.insert(0,self.node_idx[pipe.node_1].point)
            pipe.points.append(self.node_idx[pipe.node_2].point)

#    def Demand_Monte_Carlo(self,juncs_to_use,monte_variables,iterations):
#    		if len(monte_variables) == 1:
#    			monte_variables = np.ones(len(juncs_to_use))*monte_variables

#		    ##	Load the files
#		    	lib_epanet = CDLL(ENModel.dll_loc)
#		    self.lib_epanet = CDLL(ENModel.dll_loc)
#            		rep_filename = self.inp_filename.lower().rstrip('.inp') + '.rep'
#		    out_filename = self.inp_filename.lower().rstrip('.inp') + '.out'
#		    lib_epanet.ENopen(self.inp_filename, rep_filename, out_filename)
#		
#		    ###Setting it up to run a hydraulic model	
#		    lib_epanet.ENopenH()
#            lib_epanet.ENsettimeparam(0, 0)#24*3600) # sim duration = 24h
#            #lib_epanet.ENsettimeparam(1, 3600/4)  # hyd timestep = 0.25h
#            #lib_epanet.ENsettimeparam(5, 3600/4)  # reporting timestep = 0.25h         
#            lib_epanet.ENinitH(0)

#		### Some variables that we will use later
#        	time, time_step = c_long(), c_long()
#        	i, num_pipes, num_nodes = c_int(), c_int(), c_int()
#        	flow, abs_velocity= c_float(), c_float()
#        	head,demand = c_float(), c_float()
#        	pipe_id = create_string_buffer(32)
#        	node_id = create_string_buffer(32)
#        	real_pipe_id = ""
#        	velocity = 0	
#        	
#        	###########
#        	###	Node Index and Node IDs are different
#        	#	The ID is the name given in the .inp file
#        	#	The Index is the number that is assigned internally in EPAnet
#        	
#        	
#        	
#        	#iterations = 10
#        	for iteration in range(0,iterations):
#        		## run a simulation
#		  	lib_epanet.ENrunH(byref(time))	
#			for index in self.node_idx.iterkeys():
#				#pipe_id = 
#				Error = lib_epanet.ENgetnodeindex(index,byref(node_id))
#				
#				#index is the node name 'human readible from the EPAnet input file
#				# this node then gets stored in the memory at the point node_id
#				
#				#print node_id.value,index
#				i = int(index)				
#				
#				####		Store the head and demand from the last time step
#				lib_epanet.ENgetnodevalue(i,10,byref(head))
#				lib_epanet.ENgetnodevalue(i,1,byref(demand))
#				print index,head.value,demand.value
#				self.node_idx[node_id.value].monte_head.append(head.value)
#				self.node_idx[node_id.value].monte_demand.append(demand.value)
#				
##				
##				####		Apply a new demand to the selected nodes
#				try:
#					monte_index = juncs_to_use.index(int(node_id.value))
#					#print 'NodeID', node_id.value, 'Index', i
#					applied_demand = np.random.exponential(monte_variables[monte_index])	
#					Error = lib_epanet.ENsetnodevalue(i,1,c_float(applied_demand))
#					#print Error
#				except:
#					True
#		self.lib_epanet = lib_epanet
				
#    def TransientInitialise(self,model_dt):
#    		lib_epanet = self.lib_epanet
#    		#print min(self.pipe_lengths)
#    		
#    		for index in self.pipe_idx.iterkeys():   
#    			pipe_id = lib_epanet.ENgetpipeid(int(index),byref(node_id))
#    			
#    			print pipe_id.value,index
#    			pipe_dx = self.pipe_idx[index].celerity * model_dt
#    			pipe_calc_points = np.arange(0,self.pipe_idx[index].length,pipe_dx)
#    			up_stream_head = 1
#    			down_sream_head = 1
#    			pipe_initial_head = pipe_calc_points 
    			#print pipe_calc_points
    			
          



if __name__ == "__main__":
	global lib_epanet
	#EPAModel = ENModel(r'anytown.inp')
	EPAModel = ENModel(r'VerySimple.inp')
	
	#Junctions_To_Apply_Monte_Carlo_To = [2,3,4]#[2,3,4,5,6,7,8,9,10,11]
	#Junction_Monte_Carlo_Variables = [.5]#[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
	#Monte_Carlo_Iterations = 22
	
	#EPAModel.Demand_Monte_Carlo(Junctions_To_Apply_Monte_Carlo_To,Junction_Monte_Carlo_Variables,Monte_Carlo_Iterations)
	
	#EPAModel.TransientInitialise(0.01)
	
	
	
	#pp.figure()
	#pp.hist(EPAModel.junction_idx['2'].monte_head,bins = np.linspace(0,100,10),normed = True,alpha= 0.5) 	
	#pp.hist(EPAModel.junction_idx['3'].monte_head,bins = np.linspace(0,100,10),normed = True,alpha= 0.5) 
	#pp.hist(EPAModel.junction_idx['4'].monte_head,bins = np.linspace(0,100,10),normed = True,alpha= 0.5) 	
	#pp.show()
