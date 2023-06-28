#### Load functions
import numpy as np
import pandas as pd

def paleoRead(paleoFile):
	#####
	# Function reads in the paleoseismic record and returns array with inter-event 
	# times. Function expects oldest inter-event time first in file.
	# Input:
	# paleoFile: path to the paleoseismic file with inter-event times listed in ascending order (oldest first)
	# Output:
	# eq_inter: array of inter-event times in ascending order
	#####
	eq_inter = np.loadtxt(paleoFile) # Load file
	return eq_inter
	
def RRead(R_file):
	#####
	# Function reads in the paleoseismic record relative strain release and returns array 
	# with strain drops. Function expects oldest earthquake first in file.
	#####
	R_vals = np.loadtxt(R_file) # Load file
	return R_vals
	
def paramRead(dataPath):
	#####
	# Function reads in the paramater file (see readme file for details) and returns
	# inputted parameters
	# Input:
	# dataPath: path to ParamFile.txt
	# Output: Dictionary based on ParamFile.txt
	#####
	paramVals = pd.read_csv(dataPath+'/ParamFile.txt')
	modelParams = str2list(paramVals['params']) # Convert params to list
	if isinstance(paramVals['R_flag'].iloc[0],np.int64): # Check if integer and convert to float
		paramVals['R_flag'] = float(paramVals['R_flag'])
	elif "[" in paramVals['R_flag'].iloc[0]: # If R needs to be calculated
		temp_list = str2list(paramVals['R_flag']) # Create a temporary list
		paramVals['R_flag'] = WC_equation(temp_list) # Convert input values to an R using Wells and Coppersmith equation
	
	#### Parameter checks
	if paramVals['iterFit'].values[0] == 'Y':
		#### Parameter check:
		if paramVals['R_flag'].values[0] == 'Free':
			print("In the parameter file iterFit is set to Y. Either an R value must be specified or magnitude, geometry, and plate rate provided.")
			sys.exit(1)	
		if modelParams[0] != 'Free':
				print("In the parameter file iterFit is set to Y. Shape parameter must be set to Free")
				sys.exit(1) 
		if modelParams[1] != 'Free':
				print("In the parameter file iterFit is set to Y. Scale parameter must be set to Free")
				sys.exit(1) 
	return {'currentYr':paramVals['currentYr'].values[0],'priorEQ':paramVals['priorEQ'].values[0],'T_Forecast':paramVals['T_Forecast'].values[0],'yrsForecast':paramVals['yrsForecast'].values[0],'z0_flag':paramVals['z0_flag'].values[0],'R_flag':paramVals['R_flag'].values[0],'incOpenInt':paramVals['incOpenInt'].values[0],'model':paramVals['model'].values[0],'modelParams':modelParams,'iterFit':paramVals['iterFit'].values[0],'dataPath':dataPath}

def str2list(strs):
	#####
	# Function converts string to list and converts numeric values to floats.
	# Input:
	# strs: string version of list
	# Output: 
	# params: params as list
	#####
	params = strs.iloc[0][1:-1].split(" ")# Convert the string to a list
	for idx in range(len(params)): 
		if not params[idx].isalpha(): # Convert strings to floats
			params[idx] = float(params[idx])
	return params

def WC_equation(vals):
	#####
	# Function to calculate the average displacement of an earthquake based on the Wells and Coppersmith (1994) equations
	# Input:
	# vals: list of parameters for calculation [Magnitude Geometry AccumulationRate]
	# Output:
	# R: Earthquake drop size
	#####
	params = {'SS':{'a':-6.32,'b':0.9},'R':{'a':-0.74,'b':0.08},'N':{'a':-4.45,'b':0.63},'All':{'a':-4.80,'b':0.69}} # Equation coefficients
	AD = 10**(params[vals[1]]['a']+params[vals[1]]['b']*vals[0]) # Calculate average displacement in meters
	R = AD/(vals[2]*.001) # Calculate R value. Plate rate converted from mm/yr to m/yr
	return R


