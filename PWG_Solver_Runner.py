import pylab as pp
import numpy as np
from PWG_EPA_tools import *


FileName = 'Networks/1Pipe.inp'
#FileName = 'Networks/hanoi2.inp'


Net = Import_EPANet_Geom(FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(100)
Net.Initialise_Linear_Kalman(0.001)
#Net.Control_Input('Trials/StochasticDemandsPWGInput.csv')
#Net.MOC_Run(1000)
##Net.MOC_Run(86350)
##Net.geom_Plot(plot_Node_Names = True)
#Net.transient_Node_Plot(['6','10','13','16','21','24','31'])


for i in range(100):
	Net.X_Vector = Net.A_Matrix.dot(Net.X_Vector)+Net.U_Vector

pp.figure()
pp.subplot(211)
pp.plot(Net.X_Vector[:Net.X_Vector.size/2])
pp.subplot(212)
pp.plot(Net.X_Vector[Net.X_Vector.size/2:])

pp.show()
