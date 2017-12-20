import numpy as np
import pylab as pp
from epanettools import epanet2 as epa
from node import *
from pipe import *
from valve import *
from pump import *
from network import *
import csv

#### This file contains a list of functions that allow interaction between EPANet and the PWG-Transient-Solver

def Import_EPANet_Geom(inp_filename):
	####
	## This function will take an EPANet.inp file and create a PWG_Trans Network object from it.
	##	After running this function all other operations with the objet should be undertaken from the object functions
	
	
	#####
	##	Checking the EPANet file is in the correct units (PWG_Trans currently only works for LPS and D-W)
	inp_file = open(inp_filename, 'r')
	for line in inp_file:
		line = line.split()
		if len(line) != 0:
			if line[0] == 'Units':
				Units = line[1]
				if Units != 'LPS':
					print 'Incorrect Units (only LPS currently implemented)'
					break
			if line[0] == 'Headloss':
				Units = line[1]
			
				if Units != 'D-W':
					print 'Incorrect Headloss (only D-W currently implemented)'
					break
	
	Nodes = []
	Pipes = []
	Valves = []
	Pumps = []
	Node_idx = {}
	Pipe_idx = {}
	Valve_idx = {}
	Pump_idx = {}
	NodeCounter = 0
	
	#####
	##	Running through the input file to input the Node Coordinates
	inp_file = open(inp_filename, 'r')
	section = None
	for line in inp_file:
      	# Ignore rows that are either blank or contain column headings
            if line.startswith(';') or line.strip() == "":
                continue
            # Detect start of section
            elif line.startswith('['):
                section = line.strip()
            else:
                vals = line.split()
                if section == '[COORDINATES]':          
                	## 	Create a new Node object with inputs (name, xPos, yPos)
                	new_Node = Node(vals[0], float(vals[1]), float(vals[2]))
                	new_Node.number = NodeCounter
			NodeCounter += 1
                	##	Adding the new node to the Nodes list and the Directory of id's
                	Nodes.append(new_Node)
                	Node_idx[new_Node.Name] = new_Node
                	
        ######
        ##	Running through the input file again to add additional information to the nodes from the junctions, reservoir and tanks sections
        inp_file = open(inp_filename, 'r')
	section = None       	
	for line in inp_file:
      	# Ignore rows that are either blank or contain column headings
            if line.startswith(';') or line.strip() == "":
                continue
            # Detect start of section
            elif line.startswith('['):
                section = line.strip()
            else: 
            	vals = line.split()	
                if section == '[JUNCTIONS]':
                	
                	try:
                		Node_idx[vals[0]].zPos = float(vals[1])
                		Node_idx[vals[0]].demand = float(vals[2])/1000.
				Node_idx[vals[0]].idx = vals[0]
		        except:
		        	#print vals[0]
#		        	new_Node = Node(vals[0], 0, 0)
#                		Nodes.append(new_Node)
#                		Node_idx[new_Node.Name] = new_Node
		        	continue
		        	
		elif section == '[RESERVOIRS]':
		 	try:	
		 		#Setting the Head at the reservoir (this won't change)
                		Node_idx[vals[0]].H_0 = float(vals[1])
                		Node_idx[vals[0]].TranH = [float(vals[1])]
                		Node_idx[vals[0]].type = 'Reservoir'
				Node_idx[vals[0]].idx = vals[0]
		        except:
		        	#print vals[0]
		        	continue
		        	
		elif section == '[TANKS]':
			try:
                		#Node_idx[vals[0]].Head = vals[1]
                		Node_idx[vals[0]].type = 'Tank'
				Node_idx[vals[0]].idx = vals[0]
		        except:
		        	#print vals[0]
		        	continue
		        	
		
				
		        	
	#####
	##	Rerunning through the file to add in the link selections
	
	inp_file = open(inp_filename, 'r')
	section = None       	
	for line in inp_file:
      	# Ignore rows that are either blank or contain column headings
            if line.startswith(';') or line.strip() == "":
                continue
            # Detect start of section
            elif line.startswith('['):
                section = line.strip()
            else: 
            	vals = line.split()	
                if section == '[PIPES]':
                	try:
		        	new_Pipe = Pipe(vals[0],Node_idx[vals[1]],Node_idx[vals[2]],vals[3],vals[4],vals[5])
		        	Pipes.append(new_Pipe)
		        	Pipe_idx[new_Pipe.Name] = new_Pipe
		        except:
		        	continue
                elif section == '[VALVES]':
                	#print vals
                	new_Valve = Valve(vals[0], Node_idx[vals[1]], Node_idx[vals[2]], vals[3], vals[4], vals[5])
                	
                    	Valves.append(new_Valve)
                    	Valve_idx[new_Valve.Name] = new_Valve
                    	
                elif section == '[PUMPS]':
                	#print vals
			new_Pump = Pump(vals[0], Node_idx[vals[1]], Node_idx[vals[2]], vals[3])
			Pumps.append(new_Pump)
                    	Pump_idx[new_Pump.Name] = new_Pump
                    	
                    	
                	 	
	####  Creating the Network Class from the selections created.
	Net = Network(inp_filename.split('/')[-1][:-4], Nodes, Pipes, Valves, Pumps, Node_idx, Pipe_idx, Valve_idx, Pump_idx)
	Net.filename = inp_filename
	Net.Friction_Units = Units
	for pipe in Net.links:
		pipe.Friction_Units = Units
	return Net


#####
##	Function to determine the EPANet Error
def err(e):
    if(e>0):
        print e, et.ENgeterror(e,25)
        exit(5)

def Import_EPANet_Results(inp_filename,Network = None):
	if Network == None:
		Network = Import_EPANet_Geom(inp_filename)
	
	ret = epa.ENopen(inp_filename,inp_filename[:-3]+'rep',inp_filename[:-3]+'out')
	
	####	Opening the hydraulics results
	err(epa.ENopenH())
	err(epa.ENinitH(0))
	
	####	Running the Hydraulics Solver
	epa.ENrunH()
	
	####	Returning the head solutions (for the first time step, the only one accessible at the moment)
	ret,no_nodes=epa.ENgetcount(epa.EN_NODECOUNT)
	for index in range(1,no_nodes):
		ret,idx=epa.ENgetnodeid(index)
		ret,H0=epa.ENgetnodevalue(index, epa.EN_HEAD )
		try:
			#print Network.node_idx[idx].Name,idx
			Network.node_idx[idx].H_0 = H0
		except:
			print 'Problem getting Head for Node:', idx
			continue
			
	####	Returning the flow solutions (for the first time step, the only one accessible at the moment)
	ret,no_links=epa.ENgetcount(epa.EN_LINKCOUNT)
	for index in range(1,no_nodes):
		ret,idx=epa.ENgetlinkid(index)
		ret,Q0=epa.ENgetlinkvalue(index, epa.EN_FLOW )
		ret,V0=epa.ENgetlinkvalue(index, epa.EN_VELOCITY )
		try:
			
			Network.link_idx[idx].Q_0 = Q0
			Network.link_idx[idx].V_0 = V0
		except:
			print 'Problem getting Flow or Velocity for link:', idx
			continue
	

def EPA_inp_writer(Filename,Title,Junctions,Reservoirs,Pipes,Nodes):
	
	### This function will write an EPANet .inp file from input data supplied
	f = open(Filename,'wb')
	writer = csv.writer(f,dialect ='excel-tab')
	
	f.write('[TITLE] \n')
	f.write(Title + '\n')
	
	f.write('\n')
	f.write('[JUNCTIONS] \n')
	f.write(';ID              	Elev        	Demand \n')
	writer.writerows(Junctions)
	
	f.write('\n')
	f.write('[RESERVOIRS] \n')
	f.write(';ID              	Head        	Pattern          \n')
	writer.writerows(Reservoirs)
	
	f.write('\n')
	f.write('[PIPES] \n')
	f.write(';ID	Node1	Node2	Length	Diameter	Roughness	MinorLoss	Status \n')
	writer.writerows(Pipes)
		
	f.write('\n')
	f.write('[COORDINATES] \n')
	f.write(';Node	X-Coord	Y-Coord \n')
	writer.writerows(Nodes)
	
	f.write('\n')
	f.write('[OPTIONS] \n')
 	writer.writerows(np.array([['Units','LPS'],['Headloss','D-W'],['Specific Gravity','1'],['Viscosity','1'],['Trials','40'],['Accuracy','0.001'] ,['Unbalanced','Continue 10'],['Pattern','1'],['Demand Multiplier','1.0'],['Emitter Exponent','0.5'],['Diffusivity','1'],['Tolerance','0.01']]))
	
	f.write('\n')
	f.write('[REPORT] \n')
	f.write('Nodes All \n')
	f.write('Links All  \n')
	
	
	f.write('\n')
	f.write('[END] \n')
	
	f.close()
	
#if __name__ == '__main__':
#	Filename = 'Dave.inp'
#	Title = 'Daves Net'
#	Junctions = np.array([[1, 0, 0.01],[2, 0, 0.01]],dtype = '|S6')
#	Reservoirs = np.array([[3, 100, ]])
#	Pipes = np.array([[1,1,2,10,100,0.0001,0,'Open'],[2,2,3,10,100,0.0001,0,'Open']])
#	Nodes = np.array([[1,0,0],[2,10,0],[3,20,0]])
#	
#	EPA_inp_writer(Filename,Title,Junctions,Reservoirs,Pipes,Nodes)
