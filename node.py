import numpy as np
#import pylab as pp

class Node:
	def __init__(self,node_Name,xPos,yPos,zPos = 0,Demand = 0.0):
				
		#self.PipeID = PipeID
		self.type = "Node"
		self.Name = node_Name
		self.xPos = float(xPos)
		self.yPos = float(yPos)
		self.zPos = float(zPos)
		
		self.locationMetric = 1
		
		self.demand = Demand
		self.Transient_Demands = None
		###	Initial Steady State Head from EPAnet to be intitalised you need to run Import_EPANet_Results()
		self.H_0 = 0.
		self.pipes = []
		
		self.pipesIn = []
		self.pipesOut = []
		self.External_Flow = False
		#self.ResultsFile = open('Trials/Hanoi/'+self.Name +'.csv','wb')
		
		

	def initialise_Pressure_Dependent_Demands(self):
		self.Lp = self.demand / (np.sqrt(2*9.81*self.H_0))
			
				
			
			
	def MOC_iter(self,epsilon,time,dt):
		if self.Transient_Demands == None:
			demand = self.demand
		else:
			demand = self.Transient_Demands[int(time/dt)]
			#print time/dt, demand
		
		self.NewH = 0.
		
		self.BC = 0.
		self.CC = 0.
		for i in self.pipesIn:
			i.BPi = i.B + epsilon*i.R[-1]*abs(i.TranQ[-2,-2])
			i.CPi = i.TranH[-2,-2] + i.TranQ[-2,-2]*(i.B - i.R[-1]*abs(i.TranQ[-2,-2])*(1-epsilon))
			self.BC += 1./i.BPi
			self.CC += i.CPi / i.BPi
		
		for i in self.pipesOut:
			i.BMi = i.B + epsilon*i.R[0]*abs(i.TranQ[-2,1])
			i.CMi =  i.TranH[-2,1] - i.TranQ[-2,1]*(i.B - i.R[0]*abs(i.TranQ[-2,1])*(1-epsilon))
			self.BC += 1./i.BMi
			self.CC += i.CMi / i.BMi
			
		self.BC = 1./self.BC
		self.CC = self.BC*self.CC

		if self.External_Flow == True:
			#print Dave
			PressureFactor = (self.TranH[-1]/self.H_0) 
			CdA = self.CdA * PressureFactor #+ np.random.normal(0,0.000001)
			tau = 1
			s = np.sign(self.CC)							# Accounting for the non-linearity of the energy loss at the valve
			m = 0.5*(tau * CdA)**2 * s * self.BC
			n = (tau * CdA)**2 * s * -self.CC
			demand = -m + s * np.sqrt(m**2 - n)

			#print demand
		
		if self.type == 'Reservoir':
			self.NewH = self.H_0 #+np.random.normal(0,0.1)
		else:
			self.NewH = self.CC - self.BC*demand
		
		
		for i in self.pipesIn:
			i.TranH[-1,-1] = self.NewH
			i.TranQ[-1,-1] = -self.NewH / i.BPi + i.CPi / i.BPi
			
		
		for i in self.pipesOut:
			i.TranH[-1,0] = self.NewH
			i.TranQ[-1,0] = self.NewH / i.BMi - i.CMi / i.BMi
			
		self.TranH.append(self.NewH)	
		#if time % 100:
			
	
		
