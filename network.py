import numpy as np
import pylab as pp
from epanettools import epanet2 as epa
import time
#import networkx as nx
import scipy.sparse as ss
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
		
	

		self.nodal_CPs = {}
		
	#####
	##	Function To Input EPANet Results as Initial Conditions
	
	def Import_EPANet_Results(self):	
		ret = epa.ENopen(self.filename,self.filename[:-3]+'rep',self.filename[:-3]+'out')
		#print ret
		####	Opening the hydraulics results
		ret = epa.ENopenH()
		#print ret
		ret = epa.ENinitH(0)
		#print ret
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
				#self.link_idx[idx].R = float((Headloss/1000.) / (np.pi*Diameter/4. * Length * V0))
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
		self.CPs = 0
		self.pipes_State_Index = np.array([])

		link_count = 0
		for i in self.links:
			try:
				i.Q_0
				
			except:
				print 'No steady-state data for link', i.Name
				break

			### Ensuring the friction calc is done
			#i.LambdaCalc()		#Usng the Friciton model as calculated from EPANet manual
			i.friction = i.FF_0 	#Using data back calculated from EPANet

			self.link_lengths.append(float(i.length))	##Length of Link
			i.dx = i.c*dt			
			self.link_dx.append(i.c*dt)					##Nodal distance for link
			i.NoCPs = int(round(i.length / i.dx)+1)			##Number of CPs for link
			self.link_length_error.append(100.*abs((i.NoCPs-1)*i.dx - i.length)/i.length)	## Max discretisation error
			i.B = i.c/(9.81*i.area)					## This is the basic calc constant 
			i.R = i.friction*i.dx*abs(i.Q_0) / (2*9.81*i.diameter*i.area**2)	## This is the basic friction cal which we can use if we assume that it is fixed throughout the calcualtions
			#i.R = i.R * i.dx
			i.TranCPsX = np.linspace(0,i.length,i.NoCPs)	##	Generating the CPs
			pipe_CP_index = np.ones(i.TranCPsX.size)*link_count

			self.CPs += (i.TranCPsX).size	##working out the number of CPs
			self.pipes_State_Index = np.hstack((self.pipes_State_Index,pipe_CP_index))	##adding in the link to the state index
			
			i.CP_Node1 = self.pipes_State_Index.size - (i.TranCPsX).size
			i.CP_Node2 = self.pipes_State_Index.size -1
			i.link_number = link_count	##Adding in the link number to the link object
			link_count +=1
		
		self.X_Vector = (np.zeros((2*self.CPs+len(self.nodes)))).T #State Vector
		self.U_Vector = (np.zeros((2*self.CPs+len(self.nodes)))).T #Control Vector  (contains the demand values at each instance if not included in the state vector (which they could be and it would remain linear)  <- this should be step 2
		#print self.pipes_State_Index,self.pipes_State_Index.size
		

		####  	Generating the Matrices Required for the prediction steps
		self.A_Matrix = (ss.dok_matrix((2*self.CPs+len(self.nodes),2*self.CPs+len(self.nodes)))) #Dynamics Evolution Matrix
		self.P_Matrix = (ss.dok_matrix((2*self.CPs+len(self.nodes),2*self.CPs+len(self.nodes)))) #State Covariance Matrix
		self.Q_Matrix = (ss.dok_matrix((2*self.CPs+len(self.nodes),2*self.CPs+len(self.nodes)))) #Uncertainty Covariance Matrix



		###	Looping again to generate the actual values in the X_Vector and A matrix
		for i in self.links:
			self.X_Vector[i.CP_Node1:i.CP_Node2+1] = np.linspace(i.node1.H_0,i.node2.H_0,i.NoCPs).reshape(1,i.NoCPs)	### The inital Head in the state vector
			self.X_Vector[self.CPs+i.CP_Node1:self.CPs+i.CP_Node2+1] = i.Q_0 #The initial Flow vectors
			
			#### The main sections along each pipe
			for CP in range(i.CP_Node1+1,i.CP_Node2):
				#print 'Central',CP,self.CPs,CP+self.CPs+1
				#Hp,Ha
				self.A_Matrix[CP,CP-1] = 0.5
				#Hp,Hb
				self.A_Matrix[CP,CP+1] = 0.5
				#Hp,Qa
				self.A_Matrix[CP,CP+self.CPs-1] = 0.5*i.B - 0.5*i.R
				#Hp,Qb
				self.A_Matrix[CP,CP+self.CPs+1] = -0.5*i.B + 0.5*i.R
				
				#Qp,Ha
				self.A_Matrix[CP+self.CPs,CP-1] = 1./(2.*i.B)
				#Qp,Hb
				self.A_Matrix[CP+self.CPs,CP+1] = -1./(2.*i.B)
				#Qp,Qa
				self.A_Matrix[CP+self.CPs,CP+self.CPs-1] = 0.5*(1-i.R/i.B)
				#Qp,Qb
				self.A_Matrix[CP+self.CPs,CP+self.CPs+1] = 0.5*(1-i.R/i.B)
		
		for i in self.nodes:
			
			###	Adding the CPs of the nodes to store the nodal head and flow values at the junctions
			try:
				i.nodal_CP = i.pipesIn[0].CP_Node2
			except:
				i.nodal_CP = i.pipesOut[0].CP_Node1
			self.nodal_CPs[i.idx] = i.nodal_CP
			#print i.type
			Bc = 0
			Qext = i.demand
			for k in i.pipesOut:
				#print 'pOut',k.B
				Bc += 1./k.B
			for k in i.pipesIn:
				#print 'pIn',k.B
				Bc += 1./k.B
			Bc = 1./Bc
			#print Bc,k.B

			####
			##	Adding in the demand to the state vectors
			#print Qext
			self.X_Vector[2*self.CPs+i.number] = Qext
			self.A_Matrix[2*self.CPs+i.number,2*self.CPs+i.number] = 1
			#print self.X_Vector
			
			if i.type == 'Reservoir':
				#print 'In',i.pipesIn,'Out',i.pipesOut
				for k in i.pipesIn:
					#print 'Reservoir pipesIn',k.CP_Node2,k.CP_Node2+self.CPs,
					#Hp,Hp
					self.A_Matrix[k.CP_Node2,k.CP_Node2] = 1
					#Qp,Hp
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2] = -1./k.B
					#Qp,Ha
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2-1] = 1./k.B
					#Qp,Qa
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2+self.CPs-1] = 1. - k.R/k.B

				for k in i.pipesOut:
					#print 'Reservoir pipesOut',k.CP_Node1,k.CP_Node1+self.CPs,1. - k.R/k.B
					#Hp,Hp	
					self.A_Matrix[k.CP_Node1,k.CP_Node1] = 1
					#Qp,Hp
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1] = 1./k.B
					#Qp,Hb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+1] = -1./k.B
					#Qp,Qb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+self.CPs+1] = 1. - k.R/k.B

			if i.type == 'Node':
				#print 'In',i.pipesIn,'Out',i.pipesOut
				for k in i.pipesIn:	##Creation of the Cc element required for both Hp and Qp
					for j in i.pipesIn:
						#print 'Node pipesIn',k.CP_Node2,j.CP_Node2,j.CP_Node2-1,j.CP_Node2+self.CPs-1
						#Hp,Ha
						self.A_Matrix[k.CP_Node2,j.CP_Node2-1] = Bc/j.B
						#Hp,Qa
						self.A_Matrix[k.CP_Node2,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)
						#Qp,Ha
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node2-1] = Bc/j.B
						#Qp,Qa
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)

					for j in i.pipesOut:
						#print 'Node pipesOut',k.CP_Node2,j.CP_Node1,j.CP_Node1+1,j.CP_Node1+self.CPs+1
						#Hp,Hb
						self.A_Matrix[k.CP_Node2,j.CP_Node1+1] = Bc/j.B
						#Hp,Qb
						self.A_Matrix[k.CP_Node2,j.CP_Node1+self.CPs+1] = -Bc*(1-j.R/j.B)
						#Qp,Hb
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node1+1] = Bc/j.B
						#Qp,Qb
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node1+self.CPs+1] = -Bc*(1-j.R/j.B)

				for k in i.pipesOut:	##Creation of the Cc element required for both Hp and Qp
					for j in i.pipesIn:
						#print 'Node pipesIn',k.CP_Node1,j.CP_Node2,j.CP_Node2-1,j.CP_Node2+self.CPs-1
						#Hp,Ha
						self.A_Matrix[k.CP_Node1,j.CP_Node2-1] = Bc/j.B
						#Hp,Qa
						self.A_Matrix[k.CP_Node1,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)
						#Qp,Ha
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node2-1] = Bc/j.B
						#Qp,Qa
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)

					for j in i.pipesOut:
						#print 'Node pipesOut',k.CP_Node1,j.CP_Node1,j.CP_Node1+1,j.CP_Node1+self.CPs+1
						#Hp,Hb
						self.A_Matrix[k.CP_Node1,j.CP_Node1+1] = Bc/j.B
						#Hp,Qb
						self.A_Matrix[k.CP_Node1,j.CP_Node1+self.CPs+1] = -Bc*(1-j.R/j.B)
						#Qp,Hb
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+1] = Bc/j.B
						#Qp,Qb
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+self.CPs+1] = -Bc*(1-j.R/j.B)

				for k in i.pipesIn:
					#print Bc,k.B,Bc * Qext / k.B
					#print 'Dave'
					#print 'Final Node pipesIn',k.CP_Node2,k.CP_Node2+self.CPs,j.CP_Node2,j.CP_Node2-1,j.CP_Node2+self.CPs-1
					#self.U_Vector[k.CP_Node2] = -Bc * Qext #Adding the control to the head node				
					#self.U_Vector[k.CP_Node2+self.CPs] = Bc * Qext /k.B#Adding the control 

					
					#Qp
					self.A_Matrix[k.CP_Node2+self.CPs,:] *= -1./k.B
#					#Qp,Ha
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2-1] += 1./k.B
##					#Qp,Qa
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2+self.CPs-1] += (1-k.R/k.B)
					
					#### Demands in the state vector
					#Hp
					self.A_Matrix[k.CP_Node2,2*self.CPs+i.number] = -Bc
					#Qa
					self.A_Matrix[k.CP_Node2+self.CPs,2*self.CPs+i.number] = Bc/k.B
					

				for k in i.pipesOut:
					#print Bc,k.B,-Bc * Qext / k.B
					#print 'Final Node pipesOut',k.CP_Node1,k.CP_Node1+self.CPs,j.CP_Node1+1,j.CP_Node1+self.CPs+1
					#self.U_Vector[k.CP_Node1] = -Bc* Qext #Adding the control to the head node				
					#self.U_Vector[k.CP_Node1+self.CPs] = -Bc * Qext / k.B #Adding the control 
					
					
					#Qp
					self.A_Matrix[k.CP_Node1+self.CPs,:] *= 1./k.B
					#Qp,Hb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+1] += -1./k.B
					#Qp,Qb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+self.CPs+1] += (1-k.R/k.B)
					
					#### Demands in the state vector
					#Hp
					self.A_Matrix[k.CP_Node1,2*self.CPs+i.number] = -Bc
					#Qp
					self.A_Matrix[k.CP_Node1+self.CPs,2*self.CPs+i.number] = -Bc/k.B
		self.TransposeA = self.A_Matrix.T

	######
	##	Generating a P and Q matrix

	def initial_Uncertainty_Head(self,UC_percent):
		## Assuming uncertainty in flow and pressure at UC_precent of the initial value
		#Adding Uncertainty to the initial heads
		self.P_Matrix[:self.CPs,:self.CPs] = np.diag(np.ones(self.CPs)*UC_percent**2)
		#self.P_Matrix[:self.CPs,:self.CPs] = np.ones(self.CPs)*UC_percent**2

	def initial_Uncertainty_Flow(self,UC_percent):
		#Adding Uncertainty to the initial flows
		self.P_Matrix[self.CPs:2*self.CPs,self.CPs:2*self.CPs] = np.diag(np.ones(self.CPs)*UC_percent**2)
		#self.P_Matrix[self.CPs:2*self.CPs,self.CPs:2*self.CPs] = np.ones(self.CPs)*UC_percent**2

	def initial_BC_Uncertainty_Head(self,UC_percent):
		## Currently this function adds a fixed uncertainty to the Reservoir BC's 
		for i in self.nodes:
			if i.type == 'Reservoir':
				for k in i.pipesOut:
					self.P_Matrix[k.CP_Node1,k.CP_Node1] = UC_percent**2
				for k in i.pipesIn:
					self.P_Matrix[k.CP_Node2,k.CP_Node2] = UC_percent**2

	def initial_BC_Uncertainty_Demand(self,value):
		for i in self.nodes:
			if i.type == 'Node':
				self.P_Matrix[2*self.CPs+i.number,2*self.CPs+i.number] = value**2
				#self.P_Matrix[2*self.CPs+i.number,2*self.CPs+i.number] = self.X_Vector[2*self.CPs+i.number]*UC_percent**2
		


	def additional_Uncertainty_Head(self,UC_percent):
		#Adding additional uncertainty to the heads
		self.Q_Matrix[:self.CPs,:self.CPs] = np.diag(np.ones(self.CPs)*UC_percent**2)

	def additional_Uncertainty_Flow(self,UC_percent):
		#Adding additional uncertainty to the flows		
		self.Q_Matrix[self.CPs:2*self.CPs,self.CPs:2*self.CPs] = np.diag(np.ones(self.CPs)*UC_percent**2)
		
	def additional_BC_Uncertainty_Head(self,UC_percent):
		## Currently this function adds a percentage uncertainty to the Reservoir BC's and a fixed uncertainty of 0.1 l/s to all the demands, except those which already have a demand in which case it is specified by the same percentage as the reservoir
		for i in self.nodes:
			if i.type == 'Reservoir':
				for k in i.pipesOut:
					self.Q_Matrix[k.CP_Node1,k.CP_Node1] = UC_percent**2
				for k in i.pipesIn:
					self.Q_Matrix[k.CP_Node2,k.CP_Node2] = UC_percent**2	

	def additional_BC_Uncertainty_Demand(self,value):
		for i in self.nodes:
			if i.type == 'Node':
				self.Q_Matrix[2*self.CPs+i.number,2*self.CPs+i.number] = value**2
	
	#####
	##	Initiating the measurement arrays:
#	def initiate_measurement_array(self,nodes,data,R_Matrix = np.zeros((len(nodes),len(nodes)))):
#		
#		Net.R_Matrix = R_Matrix
		
	
	
	######
	##	Iterating the arrays (simple)
	##
	##	!!!!  DEPRECIATED  !!!!
	
	def kalman_iteration(self,iterations):
		self.times = np.arange(0,iterations*self.dt,self.dt)
		
		
		self.Heads = self.X_Vector[self.nodal_CPs.values()]	#Original Data
		self.P = self.P_Matrix.diagonal()[self.nodal_CPs.values()]
		self.Demands = self.X_Vector[2*self.CPs:]
		
		for t in self.times:
			### A test transient generation function, the reservoir quickly changes head
			if t == 0.1:
				### At t= 0.1 we apply a control input to add 3l/s to a demand with 1% sigma uncertainty 
				self.U_Vector[2*self.CPs+4] = 0.003
				self.Q_Matrix[2*self.CPs+4,2*self.CPs+4] = 0.003 * 0.01**2
			if t > 0.1:
				self.U_Vector[2*self.CPs+4] = 0.0
				self.Q_Matrix[2*self.CPs+4,2*self.CPs+4] = 0.
			
			### Propagating the state and unceratainty
			self.X_Vector = self.A_Matrix.dot(self.X_Vector)+self.U_Vector
			self.P_Matrix = self.A_Matrix.dot(self.P_Matrix).dot(self.TransposeA) + self.Q_Matrix

			### Storing the head and uncertainty at the nodal position
			self.Heads = np.vstack((self.Heads,self.X_Vector[self.nodal_CPs.values()]))
			self.P = np.vstack((self.P,self.P_Matrix.diagonal()[self.nodal_CPs.values()]))
			self.Demands = np.vstack((self.Demands,self.X_Vector[2*self.CPs:]))


			
	######
	##	Plotting the time histories of the nodal pressures
	def node_Pressure_Plot(self,nodes,plot_uncertainty = 0,plot_all = 0):
		pp.figure()
		if plot_all == 0:			#Only plotting the specific nodes		
			for node in nodes:
				Index = self.nodal_CPs.keys().index(node)
				pp.plot(self.times,self.Heads[:,Index])
		
			if plot_uncertainty == 1:
				for node in nodes:
					Index = self.nodal_CPs.keys().index(node)
					pp.fill_between(self.times,self.Heads[:,Index]+np.sqrt(self.P[:,Index]),self.Heads[:,Index]-np.sqrt(self.P[:,Index]),alpha = 0.25)
		else:					#Plotting all the nodes
			for Index in range(0,len(self.nodal_CPs)):
				pp.plot(self.times,self.Heads[:,Index])
			if plot_uncertainty == 1:
				for Index in range(0,len(self.nodal_CPs)):
					pp.fill_between(self.times,self.Heads[:,Index]+np.sqrt(self.P[:,Index]),self.Heads[:,Index]-np.sqrt(self.P[:,Index]),alpha = 0.25)



		pp.ylabel('Head (m)')
		pp.xlabel('Time (s)')
		pp.show()
	
	######
	##	Plotting the demands
	def demand_Plot(self):
		pp.figure()
		for i in self.nodes:
			pp.plot(self.times,self.Demands[:,i.number],label = i.Name)
			
		pp.legend()
		pp.ylabel('Demand (m3/s)')
		pp.xlabel('Time (s)')
		pp.show()
			

	######
	##	Plotting the network configuration
	def geom_Plot(self,plot_Node_Names = 0):
		pp.figure()
		pp.title('Network Geometry')
		for i in self.nodes:
			if i.type == 'Node':
				symbol = 'o'
				size = 20
			elif i.type == 'Reservoir':
				symbol = 's'
				size = 40
			elif i.type == 'Tank':
				symbol = (5,2)
				size = 40	
				
			pp.scatter([i.xPos],[i.yPos],marker = symbol,s = size,c='k')
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
		
