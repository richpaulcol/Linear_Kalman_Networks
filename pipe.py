import numpy as np
#import pylab as pp

class Pipe:
	#ID Node1 Node2 Length Diameter Roughness MinorLoss Status
	def __init__(self,pipe_Name,Node1,Node2,Length,Diameter,Roughness,wallThickness = None,modulus = None):
		self.Name = pipe_Name
		self.node1 = Node1
		self.node2 = Node2
		
		self.node1.pipesOut.append(self)	### All pipes start as a pipeOut from a node 
		self.node2.pipesIn.append(self)		### and end at a pipeIn
	
		self.diameter = float(Diameter)	/1000.		#[m]
		self.roughness = float(Roughness) 		
		self.wallThickness = wallThickness		#[m]
		self.modulus = modulus				#[GPa]
		
		k = 2.05e9
		
		self.area = self.diameter**2 *np.pi/4.
		#self.thickness = thickness
		
		self.PE = []
		self.KE = []
		
		self.Wall_Friction = []
		self.FrictionList = []
		self.ReList = []
	
		self.x1 = self.node1.xPos
		self.x2 = self.node2.xPos
		
		self.y1 = self.node1.yPos
		self.y2 = self.node2.yPos
		
		self.z1 = self.node1.zPos
		self.z2 = self.node2.zPos

		self.FF_0 = None
		
		if Length == '':
			self.length = np.sqrt((self.x2-self.x1)**2+(self.y2-self.y1)**2+(self.z2-self.z1)**2)
		else:
			self.length = float(Length)
			
	def Reynolds(self):
		#if self.TranQ 
		self.V = self.TranQ[-1,:] / self.area
		self.Re = self.diameter * abs(self.V) / (1.004*10.**-6)
		
		#print self.Re
		
	def Friction(self):  ##Not sure if should use log or log10 EPANet looks like it uses log, although equations suggest log10
		#self.friction = np.ones(self.Re.size)
#		if self.FF_0 != None:
		if self.Friction_Units == 'D-W':
		
		
			self.friction = 0.25 / (np.log10((self.roughness/1000.) / (3.7*self.diameter) + 5.74/(self.Re**0.9))**2) 
		
		elif self.Friction_Units == 'H-W':
			self.friction = 133.89 / (abs(self.V)**(4./27.) * self.roughness**(50./27.) * self.diameter**(1./16.)+10**-10)
		
		
		self.frictionLam = 64./(self.Re+1)
		self.friction[self.Re<2000] = self.frictionLam[self.Re<2000]
#		else:
		#self.friction = np.ones(self.Re.size)*self.FF_0
		####  Need to implement the interpolation across the transition region
		#self.friction[self.Re>=2000 and self.Re<4000] = self.frictionLam[self.Re>=2000 and self.Re<4000] 
	
	def LambdaCalc(self):
		rho = 998
		mu = 1e-3
		k = self.roughness
		d = self.diameter
		A = np.pi*self.diameter**2 / 4
		Re = rho * (self.Q_0/A) *self.diameter / mu
		
		if Re<2000:
			f = 64 / Re
		elif Re>4000:
			f = 0.25 / (np.log10(k / (3.7*d) + 5.74 / (Re**0.9))**2)
		else:
			R = Re/2000.
			Y2 = k / (3.7*d) + 5.74 / (Re**0.9)
			Y3 = -0.86859 *np.log(k / (3.7*d) + 5.74 / (4000**0.9))

			FA = Y3**-2
			FB = FA*(2-0.00514215 / (Y2*Y3))

			X1 = 7*FA - FB
			X2 = 0.128 - 17*FA + 2.5*FB
			X3 = -0.128 + 13*FA - 2*FB
			X4 = R*(0.032 - 3*FA + 0.5*FB)

			f = X1 + R*(X2 + R*(X3 + X4))
	
		self.friction = abs(f)
	
			
	def Wavespeed(self):
		if wallThickness == None:
			print 'Need Wall Thickness for pipe:', self.Name
		elif modulus == None:
			print 'Need Youngs Modulus for pipe:', self.Name
		else:
			self.c = np.sqrt(1./(1000.*(1./k + self.diameter/(self.wallThickness*self.modulus))))
			
	def MOC_iter(self,epsilon,time):
		self.NewH = np.zeros(self.NoTNodes)
		self.NewQ = np.zeros(self.NoTNodes)
		
		self.Reynolds()
		self.Friction()
		
		
		self.R = (self.friction*self.dx / (2*9.81 * self.diameter *self.area**2))
		
		CpList = self.TranH[-1,:-2] + self.TranQ[-1,:-2]*(self.B- self.R[:-2] *abs(self.TranQ[-1,:-2])*(1-epsilon))
		CmList = self.TranH[-1,2:] - self.TranQ[-1,2:]*(self.B - self.R[2:] *abs(self.TranQ[-1,2:])*(1-epsilon))
		
		BpList = self.B + epsilon*self.R[:-2]*abs(self.TranQ[-1,:-2])
		BmList = self.B + epsilon*self.R[2:]*abs(self.TranQ[-1,2:])
		
		self.NewQ[1:-1] = (CpList - CmList)/(BpList + BmList)
		self.NewH[1:-1] = (CmList*BpList + BmList*CpList)/(BpList+BmList)		
		
		"""This is currently set up to not save the pressure and flow"""
		
		#self.TranQ = np.vstack((self.TranQ,self.NewQ))
		#self.TranH = np.vstack((self.TranH,self.NewH))
		
		self.TranQ[0,:] = self.TranQ[1,:]
		self.TranQ[1,:]  = self.NewQ
		self.TranH[0,:] = self.TranH[1,:]
		self.TranH[1,:]  = self.NewH
		
		
		###  Calculating the Wall Friction at the pipe mid-length
		Mid = self.NoTNodes/2
		self.ReList.append(self.Re[Mid])
		self.FrictionList.append(self.friction[Mid])
		
		self.Wall_Friction.append(self.friction[Mid]*abs(self.V[Mid])*self.V[Mid]*998./8.)
		
		
		
		
		###	Calculating the Energy at the current Time
		
		self.PE.append(998. *self.area*9.81**2 / (2*self.c**2) * self.dx * np.sum((self.TranH[1,:] - self.TranH[0,:])**2 ))
		self.KE.append(998. * self.area  / 2 * self.dx*( np.sum((self.TranQ[1,:]/self.area)**2,axis = 0)))
		
		
