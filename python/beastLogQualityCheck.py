#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
QUALITY CHECK OF BEAST LOG FILES
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2021-08-04' 
__last_update__ = '2021-08-04'

# Import libraries
import os
import sys
import re
import pandas as pd
import numpy as np
import math
from statistics import mean

## Import custom functions
from handleDirectories import checkDirectory
from geneticAndSpatialEstimates import renameColumns, getBastaRegionDict

## Import regions and log files columns dictionaries
from beastLogDictionaries import * 




def flatlining(fileName, directory, nDemes, beastModel, bssvs=True, burnIn = 10):
	"""
	Get the number of continuous parameters that are stuck over more than "nIterations" 
	consecutive iterations outside the burnin phase
		- fileName (str): name of the beast log file  
		- directory (str): directory from which to read input files
		- nDemes (str): simulation framework (3demes or 7demes)
		- beastModel (str): type of beast model (dta, dta_adjustedBF, mascot_v8, basta)
		- nIterations (int): the minimum number of consecutive sampled iterations along which
		  a parameter is stuck
		- burnIn (int): percentage of the chain length corresponding to the burnin phase  

	It returns a dataframe with the number of parameters that are stuck outside the 
	burn in phase for more than 100 consecutive sampled iterations 
	Parameters:
		- nSim: number of the simulation
		- nSeq: sample size
		- protocol: sampling protocol
		- model: beast model
		- matrix: simulation framework (3demes, 7demes/mig1)

	"""

	##############################################
	# Load MCMC chain and rename columns
	##############################################
	directory = checkDirectory(directory)

	# Log file
	log = pd.read_csv(directory + fileName, sep = "\t", comment = "#")

	# Remove burnIn
	chainLength = log.shape[0]
	burnIn = round(chainLength * burnIn / 100)
	log = log.iloc[burnIn:]

	# Remove spaces and rename columns
	log.columns = log.columns.str.replace(' ', '')

	# Get region dictionary for Basta chains
	bastaRegionDict = None
	if 'basta' in beastModel.lower():
		bastaRegionDict = getBastaRegionDict(fileName, directory)

	# Rename columns
	log = renameColumns(log, beastModel, fileName, bastaRegionDict)

	##############################################
	# Create a dataframes of lags 0 to 10 for each variable
	# and identify flat lines
	##############################################
	cols = [x for x in log.columns if "rates_" in x]
	nR = len(set([re.sub("^.*_", "", x) for x in cols]))

	if nR == 1:
		regions = cols[0].split("_")
		regions.remove("rates")
		nR = len(regions)

	posterior_flat_lines = dict()

	for col in cols:

		if bssvs:
			# Indicator variables are coded in different ways
			if any(["indicators" in x for x in log.columns]):
				col_ind = col.replace("rates", "indicators")
				col_data = log.loc[(log[col_ind] == 1), col]
				pk = mean(log[col_ind] == 1)
			else:
				col_data = log.loc[(log[col] > 0), col]
				pk = mean(log[col] > 0)

			col_data.reset_index(drop = True, inplace = True)

			# Bayes factor
			qk = (math.log(2) + nR - 1) / (nR*(nR-1))
			prior_odds = qk / (1-qk)
			if pk == 1:
				pk = 1 - (1/log.shape[0])
			posterior_odds = pk / (1-pk)
			bf = posterior_odds / prior_odds

		else:
			col_data = log[col].reset_index(drop=True) 
			bf = 3


		if bf >= 3:
			lagged_data = pd.concat([col_data, col_data.shift(), col_data.shift(2), col_data.shift(3),
									col_data.shift(4), col_data.shift(5), col_data.shift(6),
									col_data.shift(7), col_data.shift(8), col_data.shift(9),
									col_data.shift(10)], axis=1).dropna()
			lagged_data.columns = ["lag_" + str(x) for x in range(11)]
			lagged_data = lagged_data.diff(axis = 1).drop("lag_0", axis = 1)

			posterior_flat_lines[col] = sum(lagged_data.apply(lambda x: all(abs(x) < 10e-6), axis = 1))

			if (posterior_flat_lines[col] >0):
				print(nDemes + " " + beastModel + " " + fileName + " " + col + " " + str(posterior_flat_lines[col]))

	stuck_mig_rates = [k for k,v in posterior_flat_lines.items() if posterior_flat_lines[k] > 0]

	##############################################
	# If flat lining is detected for one migration
	# rate, the chain is stored in a separate text
	# file
	##############################################
	if len(stuck_mig_rates) > 0:
		temp_file = "/pasteur/sonic/homes/maylayan/MMMI_Rage/2.Figures/" + nDemes + "/flatline_" + beastModel + ".txt"

		if os.path.isfile(temp_file):
			f = open(temp_file, "a")
		else:
			f = open(temp_file, "w")
			f.write("tracer ")
			#f.write("nDemes simulation beastModel fileName\n")
		f.write("HKY_" + nDemes + "/simulation" + re.sub("^sim|_.*$", "", fileName) + "/" + beastModel + "/" + fileName + " ")
		f.close()

	return(0)








def getStuckParameters(fileName, directory, simFramework, beastModel, nIterations = 50, burnIn = 10):
	"""
	Get the number of continuous parameters that are stuck over more than "nIterations" 
	consecutive iterations outside the burnin phase
		- fileName (str): name of the beast log file  
		- directory (str): directory from which to read input files
		- simFramework (str): simulation framework (3demes or 7demes)
		- beastModel (str): type of beast model (dta, dta_adjustedBF, mascot_v8, basta)
		- nIterations (int): the minimum number of consecutive sampled iterations along which
		  a parameter is stuck
		- burnIn (int): percentage of the chain length corresponding to the burnin phase  

	It returns a dataframe with the number of parameters that are stuck outside the 
	burn in phase for more than 100 consecutive sampled iterations 
	Parameters:
		- nSim: number of the simulation
		- nSeq: sample size
		- protocol: sampling protocol
		- model: beast model
		- matrix: simulation framework (3demes, 7demes/mig1)

	"""

	##############################################
	# Load and rename the dataframe
	##############################################
	# Directory
	directory = checkDirectory(directory)

	# Log file
	log = pd.read_csv(directory + fileName, sep = "\t", comment = "#")
	
	# Remove burnIn
	chainLength = log.shape[0]
	burnIn = round(chainLength * burnIn / 100)
	log = log.iloc[burnIn:]

	# Remove spaces and rename columns
	log.columns = log.columns.str.replace(' ', '')

	# Get region dictionary for Basta chains
	bastaRegionDict = None
	if 'basta' in beastModel.lower():
		bastaRegionDict = getBastaRegionDict(fileName, directory)

	# Rename columns
	log = renameColumns(log, beastModel, fileName, bastaRegionDict)

	##############################################
	# 1st order lag on continuous variables
	##############################################
	# Discard discrete and character variables 
	colToDrop = ['rootLocation', 'sumNonZeroRates', 'allTransitions'] + \
		[c for c in log.columns if "indicators" in c or "nMigration" in c] 
	log.drop(columns=[x for x in colToDrop if x in log.columns], inplace = True)

	# Round to 6 digits
	log = log.round(6)

	# Select columns with 1st order lag at 0
	log = log.loc[:,(log.diff()==0).any()]

	# Get index where >2 parameters have a lag equal to 0
	linesWithZeroes = log.diff()[(log.diff() == 0).any(axis=1)]
	nZeroesPerLine = linesWithZeroes.apply(lambda x: sum(x == 0), axis =1)
	index = list(nZeroesPerLine[nZeroesPerLine >1].index)
	nZeroesPerLine = list(nZeroesPerLine)

	##############################################
	# Max number of consecutive zeroes
	##############################################
	# Initialization
	count = 0
	counts = []
	

	# List of identical consecutive values
	if len(index) >= nIterations:
		i_prev = index[0]
		i_curr = 0
		v_prev = nZeroesPerLine[0]
		v_curr = 0

		# Check for stuck parts
		for l in range(1,len(index)):
			i_curr = index[l]
			v_curr = nZeroesPerLine[l]
			if i_curr == i_prev + 1 and v_curr == v_prev:
				count += 1
			else :
				if count != 0:
					counts.append(count)
				count = 0
			i_prev = index[l]
			v_prev = nZeroesPerLine[l]

	##############################################
	# Output
	##############################################
	# Prepare output dataframe
	nSim = int(re.match(r'.*sim(\d*)_', fileName).group(1))
	nSeq = int(re.match(r'.*_(\d*)\.log.*', fileName).group(1))
	protocol = re.match(r'.*sim\d*_(.*)_\d*\.log.*', fileName).group(1)

	results = pd.DataFrame({
		'nSim': [nSim],
		'nSeq': [nSeq],
		'protocol': [protocol],
		'matrix': [simFramework],
		'model': [beastModel]
		})

	# Fill dataframe with the longest part of the chain that is stuck
	if len(counts) > 0 and max(counts) > nIterations:
		results["nStuckChains"] = max(counts)
		print(fileName)
	else:
		results["nStuckChain"] = 0


	return(results)










def consecutiveZeroes(col):
	"""
	Get max number of consecutive 0 in a pandas Series object
		- col: column of pandas DataFrama or pandas/numpy Series object

	It returns an integer or np.nan 
	"""
	colDiff=col.diff()
	meanUnstuck = np.mean(col[colDiff > 1e-10])

	# Initialization
	col = list(col)
	colDiff = list(colDiff)
	count = 0
	counts = []
	z_prev = col[1]
	z_curr = 0

	if (len(col) < 2):
		return(np.nan)

	for i in range(2,len(colDiff)):
		z_curr = colDiff[i]
		if z_prev < 1e-10 and z_curr < 1e-10 and col[i] < meanUnstuck*0.001:
			count += 1
		else :
			if count !=0:
				counts.append(count)
			count = 0
		z_prev = colDiff[i]

	if len(counts) > 0:
		return(max(counts))
	else:
		return(np.nan)
