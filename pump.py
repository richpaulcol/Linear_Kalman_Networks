import numpy as np
import pylab as pp

class Pump(object):
	"""Represents a pump in an EPANET model."""
	def __init__(self, pump_name, node_1, node_2, parameter,diameter = 100,Pump_Location = None,wallThickness = None,modulus = None,roughness = 10**-6):
		"""Create a new pump.
		
		Instance variables:
		pump_name -- name/ID of pump
		node_1 -- name/ID of one of two terminating nodes
		node_2 -- name/ID of other terminating node
		setting -- pump setting
		fixed_setting -- pump fixed setting
		"""
		self.Name = pump_name
		self.node1 = node_1
		self.node2 = node_2
		self.parameter = parameter
		
		## Adding this pipe to the nodal lists
		self.node1.pipesOut.append(self)
		self.node2.pipesIn.append(self)
		self.diameter = float(diameter)
		self.area = np.pi * self.diameter*2/ 4.
		self.roughness = roughness
		
		
		##	Collecting the spatial locations of the pipe ends
		self.x1 = self.node1.xPos
		self.x2 = self.node2.xPos
		
		self.y1 = self.node1.yPos
		self.y2 = self.node2.yPos
		
		self.z1 = self.node1.zPos
		self.z2 = self.node2.zPos
		
		self.pump_Location = Pump_Location
		
		self.wallThickness = wallThickness
		self.modulus = modulus
		k = 2.05e9
		
		self.length = np.sqrt((self.x2-self.x1)**2+(self.y2-self.y1)**2+(self.z2-self.z1)**2)
		
		
	def Wavespeed(self):
		if wallThickness == None:
			print 'Need Wall Thickness for pump:', self.Name
			
		elif modulus == None:
			print 'Need Youngs Modulus for pump:', self.Name
			
		else:
			self.c = np.sqrt(1./(1000.*(1./k + self.diameter/(self.wallThickness*self.modulus))))
			
	def Reynolds(self):
		self.V = self.TranQ[-1,:] / self.area
		self.Re = self.diameter * abs(self.V) / 10.**-6
		
		#print self.Re
		
	def Friction(self):  ##Not sure if should use log or log10 EPANet looks like it uses log, although equations suggest log10
		#self.friction = np.ones(self.Re.size)
		self.friction = 0.25 / (np.log10(self.roughness / (3.7*self.diameter) + 5.74/(self.Re**0.9))**2) 
		self.frictionLam = 64./(self.Re+1)
		self.friction[self.Re<2000] = self.frictionLam[self.Re<2000]
