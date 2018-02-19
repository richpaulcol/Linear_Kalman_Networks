import pylab as pp
import numpy as np
from PWG_EPA_tools import *
import time
import Transient_Iterator as TI
StartTime = time.time()
import filterpy.kalman as km



Iterations = 2000

FileName = 'Networks/1Pipe.inp'
FileName = 'Networks/hanoi2.inp'
#FileName = 'Networks/1PipeReversed.inp'
FileName = 'Networks/5PipesUnseen.inp'

Net = Import_EPANet_Geom(FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(300)
Net.Initialise_Linear_Kalman(0.05)

print 'Length Error',Net.link_length_error



Net.initial_Uncertainty_Head(0.002)  # in m
Net.initial_Uncertainty_Flow(0.001)  # in m3/s

#
Net.initial_BC_Uncertainty_Head(0.001) # % value of the initial value
Net.initial_BC_Uncertainty_Demand(0.001) # in m3/s


#Net.additional_Uncertainty_Head(0.0001)  # in m
#Net.additional_Uncertainty_Flow(0.0001)  # in m3/s

Net.additional_BC_Uncertainty_Head(0.1) # % value of the initial value
Net.additional_BC_Uncertainty_Demand(0.0001) # in m3/s

#Net.Q_Matrix = ss.dok_matrix(np.ones(Net.Q_Matrix.shape)/1000.)



fig1,axes = pp.subplots(3,4)
axes[0,0].imshow(np.sqrt(Net.P_Matrix[:Net.CPs,:Net.CPs].todense()),vmin=0)
axes[0,2].imshow(np.sqrt(Net.Q_Matrix[:Net.CPs,:Net.CPs].todense()),vmin=0)

axes[1,0].imshow(np.sqrt(Net.P_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs].todense()),vmin=0)
axes[1,2].imshow(np.sqrt(Net.Q_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs].todense()),vmin=0)

axes[2,0].imshow(np.sqrt(Net.P_Matrix[2*Net.CPs:,2*Net.CPs:].todense()),vmin=0)
axes[2,2].imshow(np.sqrt(Net.Q_Matrix[2*Net.CPs:,2*Net.CPs:].todense()),vmin=0)



#####  Running the solution
#Net.times,Net.X_Vector,Net.P_Matrix,Net.Heads,Net.Demands,Net.P = TI.kalman_iteration(Net.X_Vector,Net.A_Matrix,Net.TransposeA,Net.P_Matrix,Net.U_Vector,Net.Q_Matrix,Net.nodal_CPs,Net.CPs,Iterations,Net.dt)
#TI.forward_Prediction(Net,Iterations)

#####	Including Measurements
Measurement_Nodes = ['1','4','6']
Net.R_Matrix = np.identity(len(Measurement_Nodes))*0.05**2
Net.H_Matrix = np.zeros((len(Measurement_Nodes),Net.X_Vector.size))
for i in range(len(Measurement_Nodes)):
	Net.H_Matrix[i,Net.nodal_CPs[Measurement_Nodes[i]]] = 1

Net.TransposeH = Net.H_Matrix.transpose()
#Net.MeasurementData = np.random.uniform(99,101,Iterations)
#Net.MeasurementData = np.sin(np.arange(Iterations)*0.01)+100
Net.MeasurementData = np.vstack((np.load('MeasureData1.npy'),np.load('MeasureData4.npy'),np.load('MeasureData6.npy')))


AllMeasuredData = np.vstack((np.load('MeasureData1.npy'),np.load('MeasureData2.npy'),np.load('MeasureData3.npy'),np.load('MeasureData4.npy'),np.load('MeasureData5.npy'),np.load('MeasureData6.npy')))



#######
## Iterating

TI.kalman_iteration(Net,Iterations)
#TI.forward_Prediction(Net,Iterations)


#####
##	Plotting

#Net.node_Pressure_Plot(['UpStream','Mid','DownStream'],plot_uncertainty = 1)
Net.node_Pressure_Plot(['1','4'],plot_uncertainty = 1,plot_all = 1)
#pp.plot(Net.times,Net.MeasurementData.T[:np.size(Net.times)],'k')
pp.plot(Net.times,AllMeasuredData.T[:np.size(Net.times)],':')
pp.ylim(90,110)
#Net.node_Pressure_Plot(['4'],plot_uncertainty = 1,plot_all = 0)
#pp.plot(Net.times,Net.MeasurementData,'k')
Net.demand_Plot()
#pp.figure()
#pp.plot(Net.times,Net.Demands)
#pp.show()

axes[0,1].imshow(Net.P_Matrix[:Net.CPs,:Net.CPs].todense(),vmin=0)
axes[0,3].imshow(Net.Q_Matrix[:Net.CPs,:Net.CPs].todense(),vmin=0)

axes[1,1].imshow(np.sqrt(Net.P_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs].todense()),vmin=0)
axes[1,3].imshow(np.sqrt(Net.Q_Matrix[Net.CPs:2*Net.CPs,Net.CPs:2*Net.CPs].todense()),vmin=0)

axes[2,1].imshow(np.sqrt(Net.P_Matrix[2*Net.CPs:,2*Net.CPs:].todense()),vmin=0)
axes[2,3].imshow(np.sqrt(Net.Q_Matrix[2*Net.CPs:,2*Net.CPs:].todense()),vmin=0)


axes[0,0].set_title('P_prior')
axes[0,1].set_title('P_post')
axes[0,2].set_title('Q_prior')
axes[0,3].set_title('Q_post')

axes[0,0].set_ylabel('Heads')
axes[1,0].set_ylabel('Flows')
axes[2,0].set_ylabel('Demands')

EndTime = time.time()
print "Run Time %0.3f (s) " % (EndTime-StartTime)




f = km.KalmanFilter(dim_x = Net.X_Vector.size, dim_z = Net.Z_Vector.size)
f.x = Net.X_Vector.T
f.F = Net.A_Matrix.todense()
f.H = Net.H_Matrix
f.P = Net.P_Matrix.todense()
f.R = Net.R_Matrix
f.Q = Net.Q_Matrix.todense()

Xs = np.zeros((Iterations,f.x.size))
for i in range(int(Iterations)):
	try:
		f.update(Net.MeasurementData[:,i])
	except:
		f.x = f.x.T
	Xs[i,:] = np.array(f.x)[:,0]

for n in Net.nodes:
	pp.plot(Xs[:,n.nodal_CP])

