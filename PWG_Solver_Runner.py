import pylab as pp
import numpy as np
from PWG_EPA_tools import *


FileName = 'Trials/hanoi2.inp'

Net = Import_EPANet_Geom(FileName)
Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()


Net.Constant_Wavespeed(500)
Net.MOC_Initialisation(0.1)
Net.Control_Input('Trials/StochasticDemandsPWGInput.csv')
Net.MOC_Run(1000)
#Net.MOC_Run(86350)
#Net.geom_Plot(plot_Node_Names = True)
Net.transient_Node_Plot(['6','10','13','16','21','24','31'])


#for Node in Net.nodes:
#	np.save(Node.Name,Node.TranH)
	
#for Node in Net.nodes:
#	pp.scatter([int(Node.Name)], [np.mean(Node.TranH)])
PE = np.zeros(9999)#999)
KE = np.zeros(9999)
for pipe in Net.pipes:
	PE += np.array(pipe.PE)
	KE += np.array(pipe.KE)


