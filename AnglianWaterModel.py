import pylab as pp
import numpy as np
from PWG_EPA_tools import *
import time
import Transient_Iterator as TI
StartTime = time.time()
import filterpy.kalman as km
import sys
import pykalman as pk


Iterations = 6000
Directory = 'SimpleBranched/'


FileName = Directory + 'SimpleBranched.inp'

Net = Import_EPANet_Geom(FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(1000)
Net.Initialise_Linear_Kalman(1./100)

print 'Length Error',Net.link_length_error

Net.initial_Uncertainty_Head(0.01)  # in m
Net.initial_Uncertainty_Flow(1e-1)  # in m3/s

#
Net.initial_BC_Uncertainty_Head(0.01) # % value of the initial value
Net.initial_BC_Uncertainty_Demand(1e-1) # in m3/s


Net.additional_Uncertainty_Head(1e-1)  # in m
Net.additional_Uncertainty_Flow(1e-4)  # in m3/s

Net.additional_BC_Uncertainty_Head(10.) # % value of the initial value
Net.additional_BC_Uncertainty_Demand(1e-3) # in m3/s

Net.Initial_X = np.copy(Net.X_Vector)
Net.Initial_P = np.copy(Net.P_Matrix.todense())




#####	Including Measurements
All_Nodes = ['1','2','3','4','5','6','7','8','9','10']
Measurement_Nodes = ['1','2','6','10']
NotMeasured = [x for x in All_Nodes if x not in Measurement_Nodes]


Net.R_Matrix = np.identity(len(Measurement_Nodes))*0.1
Net.H_Matrix = np.zeros((len(Measurement_Nodes),Net.X_Vector.size))

for i in range(len(Measurement_Nodes)):
	Net.H_Matrix[i,Net.nodal_CPs[Measurement_Nodes[i]]] = 1

np.random.seed(0)
MeasurementNoise = np.random.normal(0,0.1,[len(All_Nodes),1000])


Net.TransposeH = Net.H_Matrix.transpose()


Measurements = {}
Measurements['1'] = Directory+'Node_Head_1.npy'
Measurements['2'] = Directory+'Node_Head_2.npy'
Measurements['3'] = Directory+'Node_Head_3.npy'
Measurements['4'] = Directory+'Node_Head_4.npy'
Measurements['5'] = Directory+'Node_Head_5.npy'
Measurements['6'] = Directory+'Node_Head_6.npy'
Measurements['7'] = Directory+'Node_Head_7.npy'
Measurements['8'] = Directory+'Node_Head_8.npy'
Measurements['9'] = Directory+'Node_Head_9.npy'
Measurements['10'] = Directory+'Node_Head_10.npy'


MeasurementData = []
for i in Measurement_Nodes:
	MeasurementData.append(np.load(Measurements[i]))
MeasurementData = np.array(MeasurementData).T

AllMeasurementData = []
for i in All_Nodes:
	AllMeasurementData.append(np.load(Measurements[i]))
AllMeasurementData = np.array(AllMeasurementData).T	

kf = km.KalmanFilter(dim_x = Net.X_Vector.size, dim_z = MeasurementData[1,:].size)
#f = km.SquareRootKalmanFilter(dim_x = Net.X_Vector.size, dim_z = Net.MeasurementData[:,i].size)

kf.x = Net.X_Vector.T
kf.F = Net.A_Matrix.todense()
kf.H = Net.H_Matrix
kf.P = Net.P_Matrix.todense()
kf.R = Net.R_Matrix
kf.Q = Net.Q_Matrix.todense()
#kf.alpha = 1.01


X_kf = np.zeros((Iterations,kf.x.size))
#######
## Iterating
#kf.x = kf.x.T
for i in range(0,Iterations):
	
	print i
	kf.predict()
	
	if kf.H.shape[1] == kf.x.shape[0]:
		kf.update(MeasurementData[i,:].T)
	else:
		print 'Dave'
		kf.x = kf.x.T
		kf.update(MeasurementData[i,:].T)
#	try:
#		kf.update(MeasurementData[i,:].T)
#	except:
#		print('It sort of didnt work')
#		kf.x = kf.x.T
	#kf.update(MeasurementData[i,:].T)
	X_kf[i,:] = np.array(kf.x)[:,0]
	
	
colors = ['r','g','b','y','c','k','r','g','b','y','c','k']
fig,axs = pp.subplots(3,1,sharex=True)
for i in range(len(Measurement_Nodes)):
	for j in Net.nodes:
		if j.Name == Measurement_Nodes[i]:
			axs[0].plot(MeasurementData[:,i],label = j.Name+' Measured',marker='x',linewidth=0,color = colors[i])
			axs[0].plot(X_kf[:,j.nodal_CP],label=j.Name+' Filtered',color = colors[i])
#axs[0].legend()
axs[0].set_ylabel(r'Measured and Modelled Head')

for i in range(len(NotMeasured)):
	for j in Net.nodes:
		if j.Name == NotMeasured[i]:
			axs[1].plot(AllMeasurementData[:,i-1,],label = j.Name+' Measured',marker='o',linewidth=0,color = colors[i])
			axs[1].plot(X_kf[:,j.nodal_CP],label=j.Name+' Filtered',color = colors[i])
#axs[1].legend()
#axs[1].plot(AllMeasuredData[2,:])
#axs[1].plot(X_kf[:,13])
axs[1].set_xlim(0,Iterations)
axs[1].set_ylabel('True and Modelled Head')

for i in range(len(Net.nodes)):
	axs[2].plot(X_kf[:,2*Net.CPs+i],label = Net.nodes[i].Name)
axs[2].legend()
axs[2].set_ylabel('Nodal Demands')
pp.show()
