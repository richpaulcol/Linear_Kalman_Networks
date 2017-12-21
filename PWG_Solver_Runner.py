import pylab as pp
import numpy as np
from PWG_EPA_tools import *


FileName = 'Networks/1Pipe.inp'
FileName = 'Networks/hanoi2.inp'
#FileName = 'Networks/1PipeReversed.inp'
#FileName = 'Networks/2Pipes.inp'

Net = Import_EPANet_Geom(FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(1000)
Net.Initialise_Linear_Kalman(0.05)

print 'Length Error',Net.link_length_error
X0 = Net.X_Vector

#print Net.X_Vector[Net.X_Vector.size/2:]
#Net.initiating_Uncertainty()
Net.initial_BC_Uncertainty_Only(0.01)
Net.kalman_iteration(100/Net.dt)

#Net.node_Pressure_Plot(['UpStream','Mid','DownStream'],plot_uncertainty = 1)
Net.node_Pressure_Plot(['1','4'],plot_uncertainty = 1,plot_all = 1)
Net.node_Pressure_Plot(['2'],plot_uncertainty = 1,plot_all = 0)
#A = Net.A_Matrix.todense()

