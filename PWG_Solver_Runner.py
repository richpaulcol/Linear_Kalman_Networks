import pylab as pp
import numpy as np
from PWG_EPA_tools import *
import time
import Transient_Iterator as TI
StartTime = time.time()




Iterations = 1000

FileName = 'Networks/1Pipe.inp'
FileName = 'Networks/hanoi2.inp'
#FileName = 'Networks/1PipeReversed.inp'
FileName = 'Networks/5Pipes.inp'

Net = Import_EPANet_Geom(FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(300)
Net.Initialise_Linear_Kalman(0.05)

print 'Length Error',Net.link_length_error



Net.initiating_initial_Uncertainty(0.01)
Net.initiating_Additional_Uncertainty(0.0)
Net.initial_BC_Uncertainty_Head(0.0001)
Net.initial_BC_Uncertainty_Demand(1.) # in m3/s


#####  Running the solution
#Net.times,Net.X_Vector,Net.P_Matrix,Net.Heads,Net.Demands,Net.P = TI.kalman_iteration(Net.X_Vector,Net.A_Matrix,Net.TransposeA,Net.P_Matrix,Net.U_Vector,Net.Q_Matrix,Net.nodal_CPs,Net.CPs,Iterations,Net.dt)
#TI.forward_Prediction(Net,Iterations)

#####	Including Measurements
Measurement_Nodes = ['2','5']
Net.R_Matrix = ss.identity(len(Measurement_Nodes))*0.005**2
Net.H_Matrix = ss.dok_matrix((len(Measurement_Nodes),Net.X_Vector.size))
for i in range(len(Measurement_Nodes)):
	Net.H_Matrix[i,Net.nodal_CPs[Measurement_Nodes[i]]] = 1

Net.TransposeH = Net.H_Matrix.transpose()
#Net.MeasurementData = np.random.uniform(99,101,Iterations)
#Net.MeasurementData = np.sin(np.arange(Iterations)*0.01)+100
Net.MeasurementData = np.vstack((np.load('MeasureData2.npy'),np.load('MeasureData5.npy')))

TI.kalman_iteration(Net,Iterations)

#Net.initiate_measurement_array(['4'],Measurement,)
#Net.node_Pressure_Plot(['UpStream','Mid','DownStream'],plot_uncertainty = 1)
Net.node_Pressure_Plot(['1','4'],plot_uncertainty = 1,plot_all = 1)
pp.plot(Net.times,Net.MeasurementData.T,'k')
pp.ylim(90,110)
#Net.node_Pressure_Plot(['4'],plot_uncertainty = 1,plot_all = 0)
#pp.plot(Net.times,Net.MeasurementData,'k')
Net.demand_Plot()
#pp.figure()
#pp.plot(Net.times,Net.Demands)
#pp.show()

EndTime = time.time()
print "Run Time %0.3f (s) " % (EndTime-StartTime)



