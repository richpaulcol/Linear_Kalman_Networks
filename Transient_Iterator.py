import numpy as np
#import pylab as pp
#from epanettools import epanet2 as epa
#import time
#import networkx as nx
import scipy.sparse as ss
""" 


Currently only Nodal Pressures are being saved



"""
######
##	Iterating the arrays (simple)
def kalman_iteration(X_Vector,A_Matrix,TransposeA,P_Matrix,U_Vector,Q_Matrix,nodal_CPs,CPs,iterations,dt):
	times = np.arange(0,iterations*dt,dt)
	
	iterations = int(iterations)
	
	node_CPs = np.array(nodal_CPs.values())
	Heads = np.zeros((iterations,node_CPs.size))	#Original Data
	P = np.zeros((iterations,node_CPs.size))
	
	Demands = np.zeros((iterations,X_Vector[2*CPs:].size))
	
	P[0,:] = P_Matrix.diagonal()[node_CPs]
	Heads[0,:] = X_Vector[node_CPs]
	Demands[0,:] = X_Vector[2*CPs:]
	
	A_Matrix = A_Matrix.tocsr()
	TransposeA = TransposeA.tocsr()
	P_Matrix = P_Matrix.tocsr()
	
	for t in range(iterations):
		### A test transient generation function, the reservoir quickly changes head
		#print t
		
		#print t,int(0.1/dt)
		if t == int(0.1/dt):
			### At t= 0.1 we apply a control input to add 3l/s to a demand with 1% sigma uncertainty 
			U_Vector[2*CPs+4] = 0.003
			Q_Matrix[2*CPs+4,2*CPs+4] = 0.003 * 0.01**2
		if t > 0.1/dt:
			U_Vector[2*CPs+4] = 0.0
			Q_Matrix[2*CPs+4,2*CPs+4] = 0.
		
		### Propagating the state and unceratainty
		X_Vector = A_Matrix.dot(X_Vector)+U_Vector
		P_Matrix = A_Matrix.dot(P_Matrix).dot(TransposeA) + Q_Matrix

		### Storing the head and uncertainty at the nodal position
#		Heads = np.vstack((Heads,X_Vector[node_CPs]))
#		P = np.vstack((P,P_Matrix.diagonal()[node_CPs]))
#		Demands = np.vstack((Demands,X_Vector[2*CPs:]))
		
		P[t,:] = P_Matrix.diagonal()[node_CPs]
		Heads[t,:] = X_Vector[node_CPs]
		Demands[t,:] = X_Vector[2*CPs:]
		
		
	return times,X_Vector,P_Matrix,Heads,Demands,P
	
	
######
##	Iterating the arrays with inference
def kalman_iteration(X_Vector,A_Matrix,TransposeA,P_Matrix,U_Vector,Q_Matrix,nodal_CPs,CPs,iterations,dt):
	times = np.arange(0,iterations*dt,dt)
	
	iterations = int(iterations)
	
	node_CPs = np.array(nodal_CPs.values())
	Heads = np.zeros((iterations,node_CPs.size))	#Original Data
	P = np.zeros((iterations,node_CPs.size))
	
	Demands = np.zeros((iterations,X_Vector[2*CPs:].size))
	
	P[0,:] = P_Matrix.diagonal()[node_CPs]
	Heads[0,:] = X_Vector[node_CPs]
	Demands[0,:] = X_Vector[2*CPs:]
	
	A_Matrix = A_Matrix.tocsr()
	TransposeA = TransposeA.tocsr()
	P_Matrix = P_Matrix.tocsr()
	
	for t in range(iterations):
		### A test transient generation function, the reservoir quickly changes head
		#print t
		
		#print t,int(0.1/dt)
		if t == int(0.1/dt):
			### At t= 0.1 we apply a control input to add 3l/s to a demand with 1% sigma uncertainty 
			U_Vector[2*CPs+4] = 0.003
			Q_Matrix[2*CPs+4,2*CPs+4] = 0.003 * 0.01**2
		if t > 0.1/dt:
			U_Vector[2*CPs+4] = 0.0
			Q_Matrix[2*CPs+4,2*CPs+4] = 0.
		
		### Propagating the state and unceratainty
		X_Vector = A_Matrix.dot(X_Vector)+U_Vector
		P_Matrix = A_Matrix.dot(P_Matrix).dot(TransposeA) + Q_Matrix

		### Storing the head and uncertainty at the nodal position
#		Heads = np.vstack((Heads,X_Vector[node_CPs]))
#		P = np.vstack((P,P_Matrix.diagonal()[node_CPs]))
#		Demands = np.vstack((Demands,X_Vector[2*CPs:]))
		
		P[t,:] = P_Matrix.diagonal()[node_CPs]
		Heads[t,:] = X_Vector[node_CPs]
		Demands[t,:] = X_Vector[2*CPs:]
		
		
	return times,X_Vector,P_Matrix,Heads,Demands,P
