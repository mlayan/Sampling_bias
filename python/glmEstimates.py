#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MODULE TO ANALYZE GLM OUTPUTS
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-02-28' 
__last_update__ = '2020-05-04'

# Import libraries
import os
import sys
import re
import math
import pandas as pd
import numpy as np

## Import custom functions
from mathFunctions import *
from handleDirectories import checkDirectory
from geneticAndSpatialEstimates import renameColumns, rootLocationAnalysis

## Import regions and log files columns dictionaries
from beastLogDictionaries import * 







def glmMarkovJumps(
	fileName, 
	directory, 
	nDemes, 
	regionDict,
	burnIn = 10
	):
	"""
	Function to retrieve genetic and spatial estimates from Beast Log Files
		- fileName (dataframe): name of the beast log file  
		- directory (str): directory from which to read input files
		- nDemes (int): number of the spatial coupling matrix
		- regionDict (dict): dictionary whose keys are integers and values are the corresponding 
		location strings 

	It returns a dataframe with all genetic and spatial parameters of interest as lines 
	and descritive statistics as columns
	Descriptive statistics:
		- mean
		- std
		- min
		- 2.5-quantile (2.5%) of the credibility interval
		- median (50%)
		- 97.5-quantile (97.5%) of the credibility interval
		- max
		- 2.5-quantile (2.5%_hpd) of the highest posterior density
		- 97.5-quantile (97.5%_hpd) of the highest posterior density
		- BF : Bayes Factor of migration rates when BSSVS is implemented
		- value: only for root location probabilities and Kullback_Leibler divergence of root locations
	Parameters:
		- parameter: name of the parameter
		- nSim: number of the simulation
		- nSeq: sample size
		- protocol: sampling protocol
		- matrix: number of the spatial coupling matrix 
	"""
	print(fileName)

	##############################################
	# Load and rename the dataframe
	##############################################
	# Directory
	directory = checkDirectory(directory)

	# Log file
	log = pd.read_csv(directory + fileName, sep = "\t", comment = "#")
	
	# Remove burnIn
	chainLength = log.shape[0]
	burnIn = math.ceil(chainLength * burnIn / 100)
	log = log.iloc[burnIn:]

	# Remove Ne.Region columns
	toDrop = [x for x in log.columns if "Ne.Region" in x]
	log.drop(columns=toDrop, inplace=True)

	# Rename columns
	log = renameColumns(log, "mascot", fileName)
	scalers = log.filter(regex=("scaler_.*"))
	log.drop(columns=[x for x in log.columns if "scaler" in x], inplace = True)

	##############################################
	# General informations
	##############################################
	# Information on the data structure
	nSim = int(re.match(r'.*sim(\d*)_', fileName).group(1))
	nSeq = int(re.match(r'.*_(\d*)\.log.*', fileName).group(1))
	protocol = re.match(r'.*sim\d*_(.*)_\d*\.log.*', fileName).group(1)
	
	# Number of regions 
	regions = [re.sub('^.*_', '', x) for x in log.columns.tolist() if 'forwards_' in x]
	nRegions = len(set(regions))

	# Symmetric or asymmetric matrix 
	if len(regions) == (nRegions * (nRegions-1)) / 2:
		symmetric = True
	elif len(regions) == nRegions * (nRegions - 1):
		symmetric = False
	else:
		raise ValueError("The number of migration rates ({0}) doesn't correspond to \
							the number of regions ({1})".format(len(regions), nRegions))

	# List of regions
	regions = list(set(regions))

	# Logging frequency
	stepInterval = log['state'].iloc[1] - log['state'].iloc[0]

	##############################################
	# Summary statistics and ESS of all parameters
	##############################################
	# Summary statistics (min, mean, median, max, 95% CI)
	all_percentiles = [x for x in range(5, 100, 5)] + [1, 2.5, 99, 97.5]
	all_percentiles = sorted([x/100 for x in all_percentiles])
	results = log.drop(columns=['state']).describe(percentiles =all_percentiles).transpose(copy = True)
	results.drop(['count'], axis = 1, inplace = True)

	# ESS of all continuous parameters 
	# Migration rates are excluded when BSSVS is implemented
	discrete_or_rates = [x for x in log.columns if x in ["state", "sumMigPredNonZero", "sumNePredNonZero"]] + \
						[x for x in log.columns if "nMigration" in x]
	
	results['ESS'] = log.drop(columns=discrete_or_rates).apply(
		lambda x: ESS(x, stepInterval)
		)

	# 95% HPD intervals 
	results['2_5_hpd'] = log.drop(columns=['state']).apply(
		lambda x: hpd(x, lower_only = True)
		) 
	results['97_5_hpd'] = log.drop(columns=['state']).apply(
		lambda x: hpd(x, upper_only = True)
		)

	# Change results index
	results['parameter'] = results.index
	results.reset_index(inplace = True, drop = True)

	# Inclusion of probabilities derived from scalers
	results_scalers = scalers.apply(lambda x: 1-mean(x==0)).to_frame().reset_index()
	results_scalers.columns = ["parameter", "inclusion_prob"]
	results_scalers.parameter = results_scalers.parameter.str.replace("scaler_", "").replace("Ne", "NeGLM_Clock")

	# Merge sum stats and inclusion probabilities
	results = results.merge(results_scalers, how="outer", on = "parameter")

	##############################################
	# Add root location
	##############################################
	rootResults = rootLocationAnalysis(
		fileName, 
		directory, 
		"mascot", 
		regions, 
		pd.DataFrame(), 
		regionDict,
		burnIn
		)

	results = pd.concat([results, rootResults], sort = False, ignore_index = True)

	##############################################
	# Add simulation number, sample size, protocol 
	# and matrix
	##############################################
	results.rename(columns={"2.5%":"2_5", "97.5%":"97_5"}, inplace = True)
	results.columns = results.columns.str.replace("%", "")
	results['nSim'] = nSim
	results['nSeq'] = nSeq
	results['protocol'] = protocol
	results['matrix'] = nDemes
	results['model'] = "glm"

	return(results)




