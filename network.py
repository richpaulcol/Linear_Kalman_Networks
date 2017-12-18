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
			
	#######
	##	Initiailise all links with their respective wavespeeds
	def Initialise_Wavespeed(self):
		for i in self.links:
			i.Wavespeed()
	
	
	######
	##	Network MOC Initialisation
	##
	def MOC_Initialisation(self,dt):
		##  Start with a check of the dx values and the error in pipelengths and whether there is any steady state data
		self.dt = dt
		self.link_lengths = []
		self.link_dx = []
		self.link_length_error = []
		for i in self.links:
			try:
				i.Q_0
			except:
				print 'No steady-state data for link', i.Name
				break
			
			i.length = float(i.length)
			self.link_lengths.append(i.length)
			i.dx = i.c*dt
			self.link_dx.append(i.dx)
			i.NoTNodes = int(i.length / i.dx)
			self.link_length_error.append(100.*abs(i.NoTNodes*i.dx - i.length)/i.length)
			
			
			### Generating the transient constants
			#i.B = np.ones(i.NoTNodes)*i.c/(9.81*i.area)
			i.B = i.c/(9.81*i.area)
			
			### Generating the Transient Nodal positions
			i.TranNodesX = np.linspace(0,i.length,i.NoTNodes)
			### Initialising the Transient Node Heads and Flows
			
			"""The nextnx two lines are removed to enable not saving the pipe pressure and flow"""
			#i.TranH = np.linspace(i.node1.H_0,i.node2.H_0,i.NoTNodes).reshape(1,i.NoTNodes)
			#i.TranQ = np.ones(i.TranNodesX.size).reshape(,i.NoTNodes)*i.Q_0
			
			i.TranH = np.linspace(i.node1.H_0,i.node2.H_0,i.NoTNodes).reshape(1,i.NoTNodes)
			i.TranH = np.vstack((i.TranH,i.TranH))
			i.TranQ = np.ones((2,i.TranNodesX.size))*i.Q_0
			
			
			### Setting up the pressure dependent demands for transient solutions
			if self.Pressure_Dependent == 1:
				for i in self.nodes:
					i.initialise_Pressure_Dependent_Demands()
					
		print 'Maximum length error =', max(self.link_length_error)
			
	#####
	##	Function to input the control file for the transient solver
	##	Currently only nodal flow control is available	
	##	Currently only works for numerical node names	
	def Control_Input(self,ControlFile):
		Control_Data = np.loadtxt(ControlFile,delimiter = ',')
		try:
			ControlNodes = Control_Data.shape[1]
		except:
			ControlNodes = 1
			Control_Data.resize(Control_Data.shape[0],1)
		for i in range(ControlNodes):
			self.node_idx[str(int(Control_Data[0,i]))].Transient_Demands = Control_Data[1:,i] / 1000. 
		
			
	
		
	def Assign_Emmiters(self,EmitterFile):
		
		Emitter_Data = np.loadtxt(EmitterFile,delimiter = ',',dtype='string')

		
		Emitter_Nodes = Emitter_Data[0,1:]
		CdAs = Emitter_Data[1,1:]
		for i in range(0,len(CdAs)):
			try:
				CdAs[i] = float(CdAs[i])
			except:
				CdAs[i] = None
		
		#print Emitter_Nodes, CdAs
		for i in range(len(Emitter_Nodes)):
			#print Emitter_Nodes[i], CdAs[i], CdAs[i]=='None'
			self.node_idx[str(int(Emitter_Nodes[i]))].External_Flow = True
			if CdAs[i] == 'None':
				CdA = self.node_idx[str(int(Emitter_Nodes[i]))].demand / np.sqrt(self.node_idx[str(int(Emitter_Nodes[i]))].H_0)
				#print CdA
				self.node_idx[str(int(Emitter_Nodes[i]))].CdA = CdA 

			else:
				self.node_idx[str(int(Emitter_Nodes[i]))].CdA = float(CdAs[i])
				

	#####
	##	Function to run the transient solution
	##		
	def MOC_Run(self,maxT):
		self.Transient_Times = np.arange(0,maxT,self.dt)
		real_time = time.time()
		for t in self.Transient_Times[1:]:
			
			for i in self.links:
				#i.Name
				i.MOC_iter(self.epsilon,t)
			for i in self.nodes:
				#i.Name
				i.MOC_iter(self.epsilon,t,self.dt)	
			if t % 100*self.dt == 0:
				print 'Time:',t,'Elapsed Time',time.time() - real_time
				real_time = time.time()
		print 'Time:',t,'Elapsed Time',time.time() - real_time	
			
	
		
	######
	##	Plotting the network configuration
	def geom_Plot(self,plot_Node_Names = 0):
		pp.figure()
		pp.title('Network Geometry')
		for i in self.nodes:
			if i.type == 'Node':
				symbol = 'o'
				size = 10
			elif i.type == 'Reservoir':
				symbol = 's'
				size = 20
			elif i.type == 'Tank':
				symbol = (5,2)
				size = 20	
				
			pp.scatter([i.xPos],[i.yPos],marker = symbol,s = size)
			if plot_Node_Names != 0:
				pp.annotate(i.Name,(i.xPos,i.yPos))
			
		for i in self.pipes:
			pp.plot([i.x1,i.x2],[i.y1,i.y2],'k')
			
		for i in self.valves:
			pp.plot([i.x1,i.x2],[i.y1,i.y2],'r')
			
		for i in self.pumps:
			pp.plot([i.x1,i.x2],[i.y1,i.y2],'g')
		pp.axis('equal')
		pp.show()
		
		
	######
	##	PLotting the initial steady-state solution
	
	##	BROKEN
	
	
	def steady_Plot(self,Node_Range = None,Link_Range = None):
		Nodes_Data = []
		Links_Data = []
		for i in self.nodes:
			Nodes_Data.append([i.xPos,i.yPos,i.H_0])
		for i in self.links:
			Links_Data.append([[i.x1,i.x2],[i.y1,i.y2],i.Q_0])
		Nodes_Data = np.array(Nodes_Data).astype('float')
		Links_Data = np.array(Links_Data)
		pp.figure()
		
		if Node_Range == None:
			Node_max = Nodes_Data[:,2].max()
			
			Node_Colours = Nodes_Data[:,2] / Node_max
		else:
			Node_Colours = (Nodes_Data[:,2] - Node_Range[0]) / (Node_Range[1] - Node_Range[0])
		if Link_Range == None:
			Links_max = Links_Data[:,2].max()
			
			Link_Colours = Links_Data[:,2].astype('float') / Links_max
		else:
			Link_Colours = (Links_Data[:,2].astype('float') - Link_Range[0]) / (Link_Range[1] - Link_Range[0])	
		
		Link_Colours = Link_Colours.astype('float')
	#	pp.scatter(Nodes_Data[:,0],Nodes_Data[:,1],c = Node_Colours)
		for i in range(Links_Data.shape[1]):
			pp.plot(Links_Data[i,0],Links_Data[i,1],color = Link_Colours[i])
		
		pp.axis('equal')
		pp.show()
		pp.tight_layout()
		
		
	####
	##	Function to plot the transient solution at the network nodes.  If no additional argument is passed then all nodes are plotted.  The additional argument should consist of a list of node_ids
	def transient_Node_Plot(self, node_ids = 0):
		fig = pp.figure()
		ax = pp.subplot(111)
		pp.title('Nodal Head')
		if node_ids == 0:
			for i in self.node_idx:
				ax.plot(self.Transient_Times,self.node_idx[i].TranH,label = i)
		else:	
			for i in node_ids:
				ax.plot(self.Transient_Times,self.node_idx[i].TranH,label = i)
		ax.legend(loc='upper center',ncol=8)
		pp.xlabel('Time (s)')
		pp.ylabel('Head (m)')
		pp.show()
		pp.tight_layout()	
				
	def best_Logger_Location(self,plot_Node_Names = 0,No_to_highlight = 0,No_Repeat_Length = 0):
		locationMetrics = []
		for i in self.nodes:
			locationMetrics.append([i,i.locationMetric])
		locationMetrics = np.array(locationMetrics)
		#print locationMetrics
		#locationMetrics.append(min(locationMetrics)/1.25)
		pp.figure()
		pp.title('Best Logger Locations')
		for i in self.nodes:
			if i.type == 'Node':
				symbol = 'o'
				size = 100
			elif i.type == 'Reservoir':
				symbol = 's'
				size = 200
			elif i.type == 'Tank':
				symbol = (5,2)
				size = 200	
				
			colour = (i.locationMetric - min(locationMetrics[:,1]))/( max(locationMetrics[:,1]) - min(locationMetrics[:,1]))
#			print colour
			pp.scatter([i.xPos],[i.yPos],marker = symbol,s = size,c = 'r',alpha = colour)
			if plot_Node_Names != 0:
				pp.annotate(i.Name,(i.xPos,i.yPos))
			
		for i in self.pipes:
			pp.plot([i.x1,i.x2],[i.y1,i.y2],'k')
			
		for i in self.valves:
			pp.plot([i.x1,i.x2],[i.y1,i.y2],'r')
			
		for i in self.pumps:
			pp.plot([i.x1,i.x2],[i.y1,i.y2],'g')
		#print locationMetrics
		Sorted = locationMetrics[locationMetrics.argsort(axis=0)[:,1],:]
		#print Sorted
		ShortestPathLengths = nx.shortest_path_length(self.Graph,weight = 'length')	
		pp.scatter([Sorted[-1][0].xPos],[Sorted[-1][0].yPos],s=200, facecolors='none', edgecolors='k',linewidth = 4)
		#print Sorted[-1][0].Name
		for k in range(1,No_to_highlight):
			#print Sorted[-k-1][0].Name
			pp.scatter([Sorted[-k-1][0].xPos],[Sorted[-k-1][0].yPos],s=200, facecolors='none', edgecolors='k',linewidth = 4)
			
			
			
			#pp.annotate(str(k+1),(Sorted[-k-1][0].xPos,Sorted[-k-1][0].yPos))
			
		#if No_to_highlight >0:
		#	self.Graph
		
			
		pp.axis('equal')
		pp.show()
