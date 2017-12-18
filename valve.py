import numpy as np
import pylab as pp

class Valve(object):
	"""Represents a valve in an EPANET model."""
	def __init__(self, valve_name, node_1, node_2, diameter, setting, fixed_setting, Valve_Location = None,wallThickness = None,modulus = None,roughness = 10**-6):
		"""Create a new valve.
		
		Instance variables:
		valve_name -- name/ID of valve
		node_1 -- name/ID of one of two terminating nodes
		node_2 -- name/ID of other terminating node
		setting -- valve setting
		fixed_setting -- valve fixed setting
		"""
		self.Name = valve_name
		self.node1 = node_1
		self.node2 = node_2
		self.setting = setting
		self.fixed_setting = fixed_setting
		self.diameter = float(diameter)/1000.
		self.area = np.pi * self.diameter*2/ 4.
		self.roughness = roughness
		
		## Adding this pipe to the nodal lists
		self.node1.pipesOut.append(self)
		self.node2.pipesIn.append(self)
		
		
		##	Collecting the spatial locations of the pipe ends
		self.x1 = self.node1.xPos
		self.x2 = self.node2.xPos
		
		self.y1 = self.node1.yPos
		self.y2 = self.node2.yPos
		
		self.z1 = self.node1.zPos
		self.z2 = self.node2.zPos
		
		self.valve_location = Valve_Location
		self.wallThickness = wallThickness
		self.modulus = modulus
		k = 2.05e9
		
		self.length = np.sqrt((self.x2-self.x1)**2+(self.y2-self.y1)**2+(self.z2-self.z1)**2)
		
	def Wavespeed(self):
		if wallThickness == None:
			print 'Need Wall Thickness for valve:', self.Name
			
		elif modulus == None:
			print 'Need Youngs Modulus for valve:', self.Name
			
		else:
			self.c = np.sqrt(1./(1000.*(1./k + self.diameter/(self.wallThickness*self.modulus))))
			
	def Reynolds(self):
		self.V = self.TranQ[-1,:] / self.area
		self.Re = self.diameter * abs(self.V) / 10.**-6
		
		#print self.Re
		
	def Friction(self):  ##Not sure if should use log or log10 EPANet looks like it uses log, although equations suggest log10
		#self.friction = np.ones(self.Re.size)
		if self.Friction_Units == 'D-W':
			self.friction = 0.25 / (np.log10((self.roughness/1000.) / (3.7*self.diameter) + 5.74/(self.Re**0.9))**2) 
			
		elif self.Friction_Units == 'H-W':
			self.friction = 133.89 / (abs(self.V)**(4./27.) * self.roughness**(50./27.) * self.diameter**(1./16.)+10**-10)
		
		self.frictionLam = 64./(self.Re+1)
		self.friction[self.Re<2000] = self.frictionLam[self.Re<2000]
