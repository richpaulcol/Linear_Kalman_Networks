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
Net.Constant_Wavespeed(100)
Net.Initialise_Linear_Kalman(0.1)

print 'Length Error',Net.link_length_error

Plot = True
if Plot == True:
	pp.figure()
	pp.subplot(211)
	pp.plot(Net.X_Vector[:Net.X_Vector.size/2])
	pp.subplot(212)
	pp.plot(Net.X_Vector[Net.X_Vector.size/2:])


print Net.X_Vector[Net.X_Vector.size/2:]
X0 = Net.X_Vector

for i in range(10):
	Net.X_Vector = Net.A_Matrix.dot(Net.X_Vector)+Net.U_Vector
Net.X_Vector[0] = 90
for j in range(10):
	for i in range(100):
		Net.X_Vector = Net.A_Matrix.dot(Net.X_Vector)+Net.U_Vector

	#print Net.X_Vector[Net.X_Vector.size/2:]

	if Plot == True:
		pp.subplot(211)
		pp.plot(Net.X_Vector[:Net.X_Vector.size/2])
		pp.subplot(212)
		pp.plot(Net.X_Vector[Net.X_Vector.size/2:])


		pp.show()

A = Net.A_Matrix.todense()
#summ = 0
#for i in range(5):
#	print A[3,i]*X0[i]
#	summ += A[3,i]*X0[i]
#	#print summ
