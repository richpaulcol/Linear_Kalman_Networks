import pylab as pp
import numpy as np
from PWG_EPA_tools import *


FileName = 'Networks/1Pipe.inp'
FileName = 'Networks/hanoi2.inp'
#FileName = 'Networks/1PipeReversed.inp'
FileName = 'Networks/2Pipes.inp'

Net = Import_EPANet_Geom(FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(100)
Net.Initialise_Linear_Kalman(0.25)


pp.figure()
pp.subplot(211)
pp.plot(Net.X_Vector[:Net.X_Vector.size/2])
pp.subplot(212)
pp.plot(Net.X_Vector[Net.X_Vector.size/2:])




for i in range(1):
	Net.X_Vector = Net.A_Matrix.dot(Net.X_Vector)+Net.U_Vector


pp.subplot(211)
pp.plot(Net.X_Vector[:Net.X_Vector.size/2])
pp.subplot(212)
pp.plot(Net.X_Vector[Net.X_Vector.size/2:])
#pp.ylim(-1,1)

pp.show()
