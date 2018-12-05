import pylab as pp
import numpy as np
from PWG_EPA_tools import *
import time
import Transient_Iterator as TI
StartTime = time.time()
import filterpy.kalman as km
import sys
import pykalman as pk


Iterations = 1000

FileName = 'Networks/1Pipe.inp'
FileName = 'Networks/hanoi2.inp'
#FileName = 'Networks/1PipeReversed.inp'
#FileName = 'Networks/Net3b.inp'

Net = Import_EPANet_Geom(FileName)
Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(500)
Net.Initialise_Linear_Kalman(0.1)

print 'Length Error',Net.link_length_error



Net.initial_Uncertainty_Head(1.)  # in m
Net.initial_Uncertainty_Flow(0.01)  # in m3/s

#
Net.initial_BC_Uncertainty_Head(1.) # % value of the initial value
Net.initial_BC_Uncertainty_Demand(0.01) # in m3/s


Net.additional_Uncertainty_Head(0.01)  # in m
Net.additional_Uncertainty_Flow(0.00001)  # in m3/s

Net.additional_BC_Uncertainty_Head(0.1) # % value of the initial value
Net.additional_BC_Uncertainty_Demand(0.00001) # in m3/s

#Net.Q_Matrix = ss.dok_matrix(np.ones(Net.Q_Matrix.shape)/1000.)
Net.P_Matrix[-1,-1] = 1.2e-8
#Net.P_Matrix[-2,-1] = 1.e-3
#Net.P_Matrix[-1,-2] = 1.e-3

Net.Initial_X = np.copy(Net.X_Vector)
Net.Initial_P = np.copy(Net.P_Matrix.todense())


try:
	pp.cholesky(Net.P_Matrix.todense())
except:
	print "Covariance Matrix not positive definite"
#	sys.exit()

#fig1,axes = pp.subplots(3,2)
#axes[0,0].imshow((Net.P_Matrix[:Net.CPs,:Net.CPs].todense()),vmin=0)
##axes[0,2].imshow((Net.Q_Matrix[:Net.CPs,:Net.CPs].todense()),vmin=0)

#axes[1,0].imshow((Net.P_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs].todense()),vmin=0)
##axes[1,2].imshow((Net.Q_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs].todense()),vmin=0)

#axes[2,0].imshow((Net.P_Matrix[2*Net.CPs:,2*Net.CPs:].todense()),vmin=0)
##axes[2,2].imshow((Net.Q_Matrix[2*Net.CPs:,2*Net.CPs:].todense()),vmin=0)



#####  Running the solution
#Net.times,Net.X_Vector,Net.P_Matrix,Net.Heads,Net.Demands,Net.P = TI.kalman_iteration(Net.X_Vector,Net.A_Matrix,Net.TransposeA,Net.P_Matrix,Net.U_Vector,Net.Q_Matrix,Net.nodal_CPs,Net.CPs,Iterations,Net.dt)
TI.forward_Prediction(Net,Iterations)

#####	Including Measurements
#Measurement_Nodes = ['1','2','4','5','6']
#Net.R_Matrix = np.identity(len(Measurement_Nodes))*0.5**2
#Net.H_Matrix = np.zeros((len(Measurement_Nodes),Net.X_Vector.size))
#for i in range(len(Measurement_Nodes)):
#	Net.H_Matrix[i,Net.nodal_CPs[Measurement_Nodes[i]]] = 1

#Net.TransposeH = Net.H_Matrix.transpose()
#Net.MeasurementData = np.random.uniform(99,101,Iterations)
#Net.MeasurementData = np.sin(np.arange(Iterations)*0.01)+100
#Net.MeasurementData = np.vstack((np.load('MeasureData1.npy'),np.load('MeasureData2.npy'),np.load('MeasureData4.npy'),np.load('MeasureData5.npy'),np.load('MeasureData6.npy')))


#AllMeasuredData = np.vstack((np.load('MeasureData1.npy'),np.load('MeasureData2.npy'),np.load('MeasureData3.npy'),np.load('MeasureData4.npy'),np.load('MeasureData5.npy'),np.load('MeasureData6.npy')))

#f = km.KalmanFilter(dim_x = Net.X_Vector.size, dim_z = Net.MeasurementData[:,i].size)
#f = km.SquareRootKalmanFilter(dim_x = Net.X_Vector.size, dim_z = Net.MeasurementData[:,i].size)

#f.x = Net.X_Vector.T
#f.F = Net.A_Matrix.todense()
#f.H = Net.H_Matrix
#f.P = Net.P_Matrix.todense()
#f.R = Net.R_Matrix
#f.Q = Net.Q_Matrix.todense()



#######
## Iterating

#TI.kalman_iteration(Net,Iterations)
#TI.forward_Prediction(Net,Iterations)


#try:
#	pp.cholesky(Net.P_Matrix.todense())
#except:
#	print "Covariance Matrix not positive definite"

######
###	Plotting

#Net.node_Pressure_Plot(['UpStream','Mid','DownStream'],plot_uncertainty = 1)
#Net.node_Pressure_Plot(['1','4'],plot_uncertainty = 1,plot_all = 1)
#pp.plot(Net.times,Net.MeasurementData.T[:np.size(Net.times)],'k')
#pp.plot(Net.times,AllMeasuredData.T[:np.size(Net.times)],':')
#pp.ylim(95,105)
#Net.node_Pressure_Plot(['4'],plot_uncertainty = 1,plot_all = 0)
#pp.plot(Net.times,Net.MeasurementData,'k')
#Net.demand_Plot()
#pp.figure()
#pp.plot(Net.times,Net.Demands)
#pp.show()

#axes[0,1].imshow(Net.P_Matrix[:Net.CPs,:Net.CPs].todense(),vmin=0)
##axes[0,3].imshow(Net.Q_Matrix[:Net.CPs,:Net.CPs].todense(),vmin=0)

#axes[1,1].imshow((Net.P_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs].todense()),vmin=0)
##axes[1,3].imshow((Net.Q_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs].todense()),vmin=0)

#axes[2,1].imshow((Net.P_Matrix[2*Net.CPs:,2*Net.CPs:].todense()),vmin=0)
##axes[2,3].imshow((Net.Q_Matrix[2*Net.CPs:,2*Net.CPs:].todense()),vmin=0)


#axes[0,0].set_title('P_prior')
#axes[0,1].set_title('P_post')
##axes[0,2].set_title('Q_prior')
##axes[0,3].set_title('Q_post')

#axes[0,0].set_ylabel('Heads')
#axes[1,0].set_ylabel('Flows')
#axes[2,0].set_ylabel('Demands')

EndTime = time.time()
print "Run Time %0.3f (s) " % (EndTime-StartTime)


Net.times = np.arange(0,Iterations*Net.dt,Net.dt)

#f = km.KalmanFilter(dim_x = Net.X_Vector.size, dim_z = Net.MeasurementData[:,i].size)
#f.x = Net.X_Vector.T
#f.F = Net.A_Matrix.todense()
#f.H = Net.H_Matrix
#f.P = Net.P_Matrix.todense()
#f.R = Net.R_Matrix
#f.Q = Net.Q_Matrix.todense()

#Xs = np.zeros((Iterations,f.x.size))
#Ps = np.zeros((Iterations,f.x.size,f.x.size))
#Xs[0,:] = np.array(f.x)
#Ps[0,:,:] = np.array(f.P)

#for i in range(1,int(Iterations)):
#	#f.x = f.x.T
#	f.predict()
#	try:
#		f.update(Net.MeasurementData[:,i])
#	except:
#		print('It sort of didnt wwork')
#		f.x = f.x.T
#	Xs[i,:] = np.array(f.x)[:,0]
#	Ps[i,:,:] = np.array(f.P)
#	
#pp.figure()

#for n in Net.nodes:
#	pp.plot(Net.times,Xs[:,n.nodal_CP])
#	pp.fill_between(Net.times,Xs[:,n.nodal_CP]+np.sqrt(Ps[:,n.nodal_CP,n.nodal_CP]),Xs[:,n.nodal_CP]-np.sqrt(Ps[:,n.nodal_CP,n.nodal_CP]),alpha = 0.25)

#pp.ylim(95,105)


#fig2,axes2 = pp.subplots(3,2)
#axes2[0,0].imshow((Ps[0,:Net.CPs,:Net.CPs]),vmin=0)
#axes2[1,0].imshow((Ps[0,Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs]),vmin=0)
#axes2[2,0].imshow((Ps[0,2*Net.CPs:,2*Net.CPs:]),vmin=0)

#axes2[0,1].imshow((Ps[-1,:Net.CPs,:Net.CPs]),vmin=0)
#axes2[1,1].imshow((Ps[-1,Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs]),vmin=0)
#axes2[2,1].imshow((Ps[-1,2*Net.CPs:,2*Net.CPs:]),vmin=0)



#####  Using pykalman
#kf = pk.KalmanFilter()
#kf.transition_matrices = Net.A_Matrix.todense()
#kf.observation_matrices = Net.H_Matrix
#kf.transition_covariance = Net.Q_Matrix.todense()
#kf.observation_covariance = Net.R_Matrix
#kf.initial_state_mean = Net.Initial_X
#kf.initial_state_covariance = Net.Initial_P

#filtered_state_estimates,filtered_state_covariances = kf.filter(Net.MeasurementData.T)
#smoothed_state_estimates,smoothed_state_covariances = kf.smooth(Net.MeasurementData.T)



#pp.figure()

#for n in Net.nodes:
#	pp.plot(Net.times,smoothed_state_estimates[:Net.times.size,n.nodal_CP])
#	pp.fill_between(Net.times,smoothed_state_estimates[:Net.times.size,n.nodal_CP]+np.sqrt(smoothed_state_covariances[:Net.times.size,n.nodal_CP,n.nodal_CP]),smoothed_state_estimates[:Net.times.size,n.nodal_CP]-np.sqrt(smoothed_state_covariances[:Net.times.size,n.nodal_CP,n.nodal_CP]),alpha = 0.25)
#	
#pp.ylim(95,105)

kf = pk.KalmanFilter()
kf.transition_matrices = Net.A_Matrix.todense()
kf.observation_matrices = Net.H_Matrix
kf.transition_covariance = Net.Q_Matrix.todense()
kf.observation_covariance = Net.R_Matrix
kf.initial_state_mean = Net.Initial_X
kf.initial_state_covariance = Net.Initial_P



#kf = TI.pykalmanIteration(Net,k_filter=True,k_smoother = False)
#pp.figure()
#pp.plot(Net.times,Net.filtered_state_estimates[:Net.times.size,2*Net.CPs:])

##pp.figure()
##pp.plot(Net.times,Net.smoothed_state_estimates[:Net.times.size,2*Net.CPs:])


#kf.em_vars=['transition_covariance']
#kf = kf.em(X=Net.MeasurementData.T, n_iter=2)
#kf.initial_state_mean = Net.Initial_X
#kf.initial_state_covariance = Net.Initial_P
#Net.filtered_state_estimates,Net.filtered_state_covariances = kf.filter(Net.MeasurementData.T)

#pp.figure()
#pp.plot(Net.times,Net.filtered_state_estimates[:Net.times.size,2*Net.CPs:])

#pp.figure()

#for n in Net.nodes:
#	pp.plot(Net.times,Net.filtered_state_estimates[:Net.times.size,n.nodal_CP])
#	pp.fill_between(Net.times,Net.filtered_state_estimates[:Net.times.size,n.nodal_CP]+np.sqrt(Net.filtered_state_covariances[:Net.times.size,n.nodal_CP,n.nodal_CP]),Net.filtered_state_estimates[:Net.times.size,n.nodal_CP]-np.sqrt(Net.filtered_state_covariances[:Net.times.size,n.nodal_CP,n.nodal_CP]),alpha = 0.25)
#	
#pp.ylim(95,105)
#pp.show()

