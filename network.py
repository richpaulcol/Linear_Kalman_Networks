import numpy as np
import pylab as pp
from epanettools import epanet2 as epa
import time
import networkx as nx
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
		
		
		self.KE = [0]
		self.PE = [0]
		
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
			print 2*9.81*(Headloss/1000.)*Diameter / (Length * V0**2)
			
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
			#i.LambdaCalc()
			i.friction = i.FF_0

			self.link_lengths.append(float(i.length))	##Length of Link
			i.dx = i.c*dt			
			self.link_dx.append(i.c*dt)					##Nodal distance for link
			i.NoTCPs = int(i.length / i.dx)			##Number of CPs for link
			self.link_length_error.append(100.*abs(i.NoTCPs*i.dx - i.length)/i.length)	## Max discretisation error
			i.B = i.c/(9.81*i.area)					## This is the basic calc constant 
			i.R = i.friction*i.dx*i.Q_0 / (2*9.81*i.diameter*i.area**2)	## This is the basic friction cal which we can use if we assume that it is fixed throughout the calcualtions
			i.TranCPsX = np.linspace(0,i.length,i.NoTCPs)	##	Generating the CPs
			pipe_CP_index = np.ones(i.TranCPsX.size)*link_count

			self.CPs += (i.TranCPsX).size	##working out the number of CPs
			self.pipes_State_Index = np.hstack((self.pipes_State_Index,pipe_CP_index))	##adding in the link to the state index
			
			i.CP_Node1 = self.pipes_State_Index.size - (i.TranCPsX).size
			i.CP_Node2 = self.pipes_State_Index.size -1
			i.link_number = link_count	##Adding in the link number to the link object
			link_count +=1
		
		self.X_Vector = (np.zeros((2*self.CPs))).T #State Vector
		self.U_Vector = (np.zeros((2*self.CPs))).T #Control Vector  (contains the demand values at each instance if not included in the state vector (which they could be and it would remain linear)  <- this should be step 2
		#print self.pipes_State_Index,self.pipes_State_Index.size
		

		####  	Generating the Matrices Required for the prediction steps
		self.A_Matrix = (ss.dok_matrix((2*self.CPs,2*self.CPs))) #Dynamics Evolution Matrix
		self.P_Matrix = (ss.dok_matrix((2*self.CPs,2*self.CPs))) #State Covariance Matrix
		self.Q_Matrix = (ss.dok_matrix((2*self.CPs,2*self.CPs))) #Uncertainty Covariance Matrix



		###	Looping again to generate the actual values in the X_Vector and A matrix
		for i in self.links:
			self.X_Vector[i.CP_Node1:i.CP_Node2+1] = np.linspace(i.node1.H_0,i.node2.H_0,i.NoTCPs).reshape(1,i.NoTCPs)	### The inital Head in the state vector
			self.X_Vector[self.CPs+i.CP_Node1:self.CPs+i.CP_Node2+1] = i.Q_0 #The initial Flow vectors
			
			#### The main sections along each pipe
			for CP in range(i.CP_Node1+1,i.CP_Node2):
				#print CP,self.CPs,CP+self.CPs+1
				#Hp,Ha
				self.A_Matrix[CP,CP-1] = 0.5
				#Hp,Hb
				self.A_Matrix[CP,CP+1] = 0.5
				#Hp,Qa
				self.A_Matrix[CP,CP+self.CPs-1] = 0.5*i.B - 0.5*i.R
				#Hp,Qb
				self.A_Matrix[CP,CP+self.CPs+1] = -0.5*i.B + 0.5*i.R
				
				#Qp,Ha
				self.A_Matrix[CP+self.CPs,CP-1] = 1./(2*i.B)
				#Qp,Hb
				self.A_Matrix[CP+self.CPs,CP+1] = -1./(2*i.B)
				#Qp,Qa
				self.A_Matrix[CP+self.CPs,CP+self.CPs-1] = 0.5*(1-i.R/i.B)
				#Qp,Qb
				self.A_Matrix[CP+self.CPs,CP+self.CPs+1] = 0.5*(1-i.R/i.B)
		
		for i in self.nodes:
			print i.type
			Bc = 0
			Qext = i.demand
			for k in i.pipesOut:
				#print 'pOut',k.B
				Bc += 1./k.B
			for k in i.pipesIn:
				#print 'pIn',k.B
				Bc += 1./k.B
			Bc = 1./Bc
			print Bc,k.B
			
			if i.type == 'Reservoir':
				for k in i.pipesIn:
					#Hp,Hp
					self.A_Matrix[k.CP_Node2,k.CP_Node2] = 1
					#Qp,Hp
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2] = -1./k.B
					#Qp,Ha
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2-1] = 1./k.B
					#Qp,Qa
					self.A_Matrix[k.CP_Node2+self.CPs,k.CP_Node2+self.CPs-1] = 1. - k.R/k.B

				for k in i.pipesOut:
					#Hp,Hp	
					self.A_Matrix[k.CP_Node1,k.CP_Node1] = 1
					#Qp,Hp
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1] = 1./k.B
					#Qp,Hb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+1] = -1./k.B
					#Qp,Qb
					self.A_Matrix[k.CP_Node1+self.CPs,k.CP_Node1+self.CPs+1] = 1. - k.R/k.B

			if i.type == 'Node':
				
				for k in i.pipesIn+i.pipesOut:	##Creation of the Cc element required for both Hp and Qp
					for j in i.pipesIn:
						#print k.CP_Node2,j.CP_Node2
						#Hp,Ha
						self.A_Matrix[k.CP_Node2,j.CP_Node2-1] = Bc/j.B
						#Hp,Qa
						self.A_Matrix[k.CP_Node2,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)
						#Qp,Ha
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node2-1] = Bc/j.B
						#Qp,Qa
						self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node2+self.CPs-1] = Bc*(1-j.R/j.B)

					for j in i.pipesOut:
						#print 'Hello'
						#Hp,Hb
						self.A_Matrix[k.CP_Node1,j.CP_Node1+1] = Bc/j.B
						#Hp,Qb
						self.A_Matrix[k.CP_Node1,j.CP_Node1+self.CPs+1] = Bc*(-1+j.R/j.B)
						#Qp,Hb
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+1] = Bc/j.B
						#Qp,Qb
						self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+self.CPs+1] = Bc*(-1+j.R/j.B)

				for k in i.pipesIn:
					print Bc,k.B,Bc * Qext / k.B
					#print 'Dave'
					self.U_Vector[k.CP_Node2] = -Bc * Qext #Adding the control to the head node				
					self.U_Vector[k.CP_Node2+self.CPs] = Bc * Qext /k.B#Adding the control 
					#Qp
					self.A_Matrix[k.CP_Node2+self.CPs,:] *= -1./k.B
					#Qp,Hp
					#self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node2] += -1./k.B
#					#Qp,Ha
					self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node2-1] += 1./k.B
##					#Qp,Qa
					self.A_Matrix[k.CP_Node2+self.CPs,j.CP_Node2+self.CPs-1] += (1-k.R/k.B)
					
					

				for k in i.pipesOut:
					print Bc,k.B,-Bc * Qext / k.B
					self.U_Vector[k.CP_Node1] = -Bc* Qext #Adding the control to the head node				
					self.U_Vector[k.CP_Node1+self.CPs] = -Bc * Qext / k.B #Adding the control 
					#Qp
					self.A_Matrix[k.CP_Node1+self.CPs,:] *= 1./k.B
					#Qp,Hp
					#self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node1] += 1./k.B
					#Qp,Hb
					self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+1] += -1./k.B
					#Qp,Qb
					self.A_Matrix[k.CP_Node1+self.CPs,j.CP_Node1+self.CPs+1] += (1-k.R/k.B)






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
		
