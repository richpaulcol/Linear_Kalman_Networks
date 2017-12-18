import numpy as np
import pylab as pp
from epanettools import epanet2 as epa
import time
import networkx as nx

""" 


Currently only Nodal Pressures are being saved



"""




def err(e):
    if(e>0):
        print e, epa.ENgeterror(e,25)
        #exit(5)


class Network(object):
	"""Represents a WDS Network."""
	def __init__(self,net_Name,nodes,pipes,valves,pumps,node_idx,pipe_idx,valve_idx,pump_idx):
		self.Name = net_Name
		self.nodes = nodes
		self.pipes = pipes
		self.valves = valves
		self.pumps = pumps
		
		#### Creating a list of link objects 
		self.links = self.pipes + self.valves + self.pumps
		
		self.node_idx = node_idx
		self.pipe_idx = pipe_idx
		self.valve_idx = valve_idx
		self.pump_idx = pump_idx
		
		#### Creating a dictionary of link objects id's
		self.link_idx = dict(self.pipe_idx.items() + self.valve_idx.items()+ self.pump_idx.items())
		
		
		#### Network Config Parameters
		self.Pressure_Dependent = 0
		self.epsilon = 0.85
		
		
		self.KE = [0]
		self.PE = [0]
		
	#####
	##	Function To Input EPANet Results as Initial Conditions
	
	def Import_EPANet_Results(self):	
		ret = epa.ENopen(self.filename,self.filename[:-3]+'rep',self.filename[:-3]+'out')
		print ret
		####	Opening the hydraulics results
		ret = epa.ENopenH()
		print ret
		ret = epa.ENinitH(0)
		print ret
		####	Running the Hydraulics Solver
		epa.ENrunH()
		
		####	Returning the head solutions (for the first time step, the only one accessible at the moment)
		##	Only need to do this for the node elements at the moment as reservoirs don't change
		ret,no_nodes=epa.ENgetcount(epa.EN_NODECOUNT)
		print 'Number of NODES in results file', no_nodes
		for index in range(1,no_nodes+1):
			ret,idx=epa.ENgetnodeid(index)
			ret,H0=epa.ENgetnodevalue(index, epa.EN_HEAD )
			try:
				#print Network.node_idx[idx].Name,idx
				if self.node_idx[idx].type == 'Node':
					self.node_idx[idx].H_0 = float(H0)
					self.node_idx[idx].TranH = [float(H0)]

			except:
				print 'Problem getting Head for Node:', idx
				continue
			
		####	Returning the flow solutions (for the first time step, the only one accessible at the moment)
		ret,no_links=epa.ENgetcount(epa.EN_LINKCOUNT)
		print 'Number of LINKS in results file',no_links
		for index in range(1,no_links+1):
			ret,idx=epa.ENgetlinkid(index)
			
			ret,Q0=epa.ENgetlinkvalue(index, epa.EN_FLOW )
			
			ret,V0=epa.ENgetlinkvalue(index, epa.EN_VELOCITY )

			ret,Headloss=epa.ENgetlinkvalue(index,epa.EN_HEADLOSS)
			ret,Length=epa.ENgetlinkvalue(index,epa.EN_LENGTH)
			ret,Diameter=epa.ENgetlinkvalue(index,epa.EN_DIAMETER)
			#print Headloss,Length,Diameter,V0
			#print 2*9.81*(Headloss/1000.)*Diameter / (Length * V0**2)
			
			try:
			
				self.link_idx[idx].Q_0 = float(Q0)/1000. #Convert to m^3/s
				self.link_idx[idx].V_0 = float(V0)
				self.link_idx[idx].FF_0 = float(2*9.81*(Headloss/1000.)*Diameter / (Length * V0**2))
			except:
				print 'Problem getting Flow or Velocity for link:', idx
				continue	
	######
	##
	##	Initialise all links with a constant wave speed
	def Constant_Wavespeed(self,wavespeed):
		for i in self.links:
			i.c = wavespeed
			
	#####
	##	
	##	Initialise the Linear Kalman Forward Propagation Matrices
	##	Linear Friction, Uncertainty in only the Flow and Head Vector
	def Initialise_Linear_Kalman(self,dt):
		self.dt = dt
		self.link_lengths = []
		self.link_dx = []
		self.link_length_error = []
		self.no_calc_points = 0

		for i in self.links:
			try:
				i.Q_0
			except:
				print 'No steady-state data for link', i.Name
				break

			self.link_lengths.append(float(i.length))	##Length of Link
			i.dx = i.c*dt			
			self.link_dx.append(i.c*dt)					##Nodal distance for link
			i.NoTNodes = int(i.length / i.dx)			##Number of nodes for link
			self.link_length_error.append(100.*abs(i.NoTNodes*i.dx - i.length)/i.length)	## Max discretisation error
			i.B = i.c/(9.81*i.area)
			i.TranNodesX = np.linspace(0,i.length,i.NoTNodes)
			self.no_calc_points += (i.TranNodesX).size
		
		self.X_Vector = (np.zeros((2*self.no_calc_points))).T
		
