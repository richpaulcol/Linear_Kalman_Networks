import numpy as np
#import pylab as pp
#from epanettools import epanet2 as epa
#import time
#import networkx as nx
import scipy.sparse as ss
from scipy.linalg import inv

import filterpy.kalman as km
import pykalman as pk


""" 


Currently only Nodal Pressures are being saved



"""
######
##	Iterating the arrays (simple)
def kalman_iteration_explicit(X_Vector,A_Matrix,TransposeA,P_Matrix,U_Vector,Q_Matrix,nodal_CPs,CPs,iterations,dt):
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
	
	
def forward_Prediction(Net,iterations):
	Net.times = np.arange(0,iterations*Net.dt,Net.dt)
	
	iterations = int(iterations)
	
	node_CPs = np.array(Net.nodal_CPs.values())
	
	Net.Heads = np.zeros((iterations,node_CPs.size))	#Original Data
	Net.P = np.zeros((iterations,node_CPs.size))
	Net.Demands = np.zeros((iterations,Net.X_Vector[2*Net.CPs:].size))
	
	Net.P[0,:] = Net.P_Matrix.diagonal()[node_CPs]
	Net.Heads[0,:] = Net.X_Vector[node_CPs]
	Net.Demands[0,:] = Net.X_Vector[2*Net.CPs:]
	
	Net.A_Matrix = Net.A_Matrix.tocsr()
	Net.TransposeA = Net.TransposeA.tocsr()
	Net.P_Matrix = Net.P_Matrix.tocsr()
	
	for t in range(iterations):
		### A test transient generation function, the reservoir quickly changes head
		#print t
		
		#print t,int(0.1/dt)
#		if t == int(0.1/Net.dt):
#			### At t= 0.1 we apply a control input to add 3l/s to a demand with 1% sigma uncertainty 
#			Net.U_Vector[2*Net.CPs+4] = 0.003
#			Net.Q_Matrix[2*Net.CPs+4,2*Net.CPs+4] = 0.003 * 0.01**2
#		if t > 0.1/Net.dt:
#			Net.U_Vector[2*Net.CPs+4] = 0.0
#			Net.Q_Matrix[2*Net.CPs+4,2*Net.CPs+4] = 0.
		
		### Propagating the state and unceratainty
		Net.X_Vector = Net.A_Matrix.dot(Net.X_Vector)+Net.U_Vector
		Net.P_Matrix = Net.A_Matrix.dot(Net.P_Matrix).dot(Net.TransposeA) + Net.Q_Matrix

		### Storing the head and uncertainty at the nodal position
#		Heads = np.vstack((Heads,X_Vector[node_CPs]))
#		P = np.vstack((P,P_Matrix.diagonal()[node_CPs]))
#		Demands = np.vstack((Demands,X_Vector[2*CPs:]))
		
		Net.P[t,:] = Net.P_Matrix.diagonal()[node_CPs]
		Net.Heads[t,:] = Net.X_Vector[node_CPs]
		Net.Demands[t,:] = Net.X_Vector[2*Net.CPs:]
		
		
	#return times,X_Vector,P_Matrix,Heads,Demands,P
	
	
	
######
##	Iterating the arrays with inference
def kalman_iteration(Net,iterations):
	Net.times = np.arange(0,iterations*Net.dt,Net.dt)
	
	iterations = int(iterations)
	
	node_CPs = np.array(Net.nodal_CPs.values())
	
	Net.Heads = np.zeros((iterations,node_CPs.size))	#Original Data
	Net.P = np.zeros((iterations,node_CPs.size))
	Net.Demands = np.zeros((iterations,Net.X_Vector[2*Net.CPs:].size))
	
	Net.P[0,:] = Net.P_Matrix.diagonal()[node_CPs]
	Net.Heads[0,:] = Net.X_Vector[node_CPs]
	Net.Demands[0,:] = Net.X_Vector[2*Net.CPs:]
	
	Net.A_Matrix = Net.A_Matrix.tocsr()
	Net.TransposeA = Net.TransposeA.tocsr()
	Net.P_Matrix = Net.P_Matrix.tocsr()
	
	for t in range(iterations):
		### A test transient generation function, the reservoir quickly changes head
#		print t
		
		#print t,int(0.1/dt)
#		if t == int(0.1/Net.dt):
#			### At t= 0.1 we apply a control input to add 3l/s to a demand with 1% sigma uncertainty 
#			Net.U_Vector[2*Net.CPs+4] = 0.003
#			Net.Q_Matrix[2*Net.CPs+4,2*Net.CPs+4] = 0.003 * 0.01**2
#		if t > 0.1/Net.dt:
#			Net.U_Vector[2*Net.CPs+4] = 0.0
#			Net.Q_Matrix[2*Net.CPs+4,2*Net.CPs+4] = 0.
		
		### Propagating the state and unceratainty
		Net.X_Vector = Net.A_Matrix.dot(Net.X_Vector)+Net.U_Vector
		Net.P_Matrix = Net.A_Matrix.dot(Net.P_Matrix).dot(Net.TransposeA) + Net.Q_Matrix
		
		### Updating the state and uncertainty using data
		
		Net.Z_Vector = Net.MeasurementData[:,t]
		
		#Net.K = Net.P_Matrix.dot(Net.TransposeH).dot( inv(Net.H_Matrix.dot(Net.P_Matrix).dot(Net.TransposeH) + Net.R_Matrix).todense() )
		Net.K = Net.P_Matrix.dot(Net.TransposeH).dot( inv(Net.H_Matrix.dot(Net.P_Matrix.todense()).dot(Net.TransposeH) + Net.R_Matrix) )

		Net.X_Vector = Net.X_Vector + Net.K.dot(Net.Z_Vector - Net.H_Matrix.dot(Net.X_Vector))

		#Net.P_Matrix = ss.csc_matrix(ss.identity(Net.X_Vector.size) - Net.K.dot(Net.H_Matrix).dot(Net.P_Matrix.todense()))
		#Net.P_Matrix = ss.csc_matrix(Net.P_Matrix - Net.K.dot((Net.H_Matrix.dot(Net.P_Matrix.todense().dot(Net.H_Matrix.T))) + Net.R_Matrix).dot(Net.K.T))

		Net.P_Matrix = ss.csc_matrix( (np.identity(Net.X_Vector.size) - Net.K.dot(Net.H_Matrix)).dot(Net.P_Matrix.todense()) )

		### Storing the head and uncertainty at the nodal position
#		Heads = np.vstack((Heads,X_Vector[node_CPs]))
#		P = np.vstack((P,P_Matrix.diagonal()[node_CPs]))
#		Demands = np.vstack((Demands,X_Vector[2*CPs:]))
		
		Net.P[t,:] = Net.P_Matrix.diagonal()[node_CPs]
		Net.Heads[t,:] = Net.X_Vector[node_CPs]
		Net.Demands[t,:] = Net.X_Vector[2*Net.CPs:]
		
def filterPyIteration(Net,iterations):
	f = km.KalmanFilter(dim_x = Net.X_Vector.size, dim_z = Net.Z_Vector.size)
	f.x = Net.X_Vector.T
	f.F = Net.A_Matrix.todense()
	f.H = Net.H_Matrix
	f.P = Net.P_Matrix.todense()
	f.R = Net.R_Matrix
	f.Q = Net.Q_Matrix.todense()
		
	Xs = np.zeros((Iterations,f.x.size))
	Ps = np.zeros((Iterations,f.x.size,f.x.size))
	Xs[0,:] = np.array(f.x)
	Ps[0,:,:] = np.array(f.P)

	for i in range(1,int(Iterations)):
		#f.x = f.x.T
		f.predict()
		try:
			f.update(Net.MeasurementData[:,i])
		except:
			print('It sort of didnt wwork')
			f.x = f.x.T
		Xs[i,:] = np.array(f.x)[:,0]
		Ps[i,:,:] = np.array(f.P)
	
	
def pykalmanIteration(Net,k_filter=True,k_smoother = False):
	kf = pk.KalmanFilter()
	kf.transition_matrices = Net.A_Matrix.todense()
	kf.observation_matrices = Net.H_Matrix
	kf.transition_covariance = Net.Q_Matrix.todense()
	kf.observation_covariance = Net.R_Matrix
	kf.initial_state_mean = Net.Initial_X
	kf.initial_state_covariance = Net.Initial_P
	
	if k_filter == True:
		Net.filtered_state_estimates,Net.filtered_state_covariances = kf.filter(Net.MeasurementData.T)
	if k_smoother == True:
		Net.smoothed_state_estimates,Net.smoothed_state_covariances = kf.smooth(Net.MeasurementData.T)
	
	return kf
