import pylab as pp
import numpy as np
from PWG_EPA_tools import *


FileName = 'Networks/hanoi2.inp'

Net = Import_EPANet_Geom(FileName)
#Net.geom_Plot(plot_Node_Names = True)
Net.Import_EPANet_Results()
Net.Constant_Wavespeed(500)
Net.Initialise_Linear_Kalman(0.1)
#Net.Control_Input('Trials/StochasticDemandsPWGInput.csv')
#Net.MOC_Run(1000)
##Net.MOC_Run(86350)
##Net.geom_Plot(plot_Node_Names = True)
#Net.transient_Node_Plot(['6','10','13','16','21','24','31'])



