#!/usr/bin/env python
# coding: utf-8

"""
COMPUTE ADJUSTED BAYES FACTORS
"""

# Import libraries
import os
import sys
import re
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

## Import custom functions
sys.path.append('python/')
from treePosteriorDistributionStatistics import *

## Set environment
directory = sys.argv[1] + "/"   # "7demes" or "3demes"
if sys.argv[1] == "3demes":
	regions = {
	'Region1':0, 
	'Region2':1, 
	'Region3':2
	}
else:
	regions = {
	'Region1':0, 
	'Region2':1, 
	'Region3':2, 
	'Region4':3, 
	'Region5':4, 
	'Region6':5, 
	'Region7':6
	}

os.chdir(directory)

# Get log file names
beastModel = sys.argv[2]	# "dta", "basta", "mascot" or "glm"
nSeq = sys.argv[3]			# "150" or "150"
logFiles = []
directories = []

for root, subDirs, files in os.walk('.'):
	if root.endswith(beastModel) and 'mig' not in root :
		for f in files:
			if 'files' in beastModel:
				if f.endswith(nSeq + '.nex'):
					logFiles.append(f)
					directories.append(root)
			else:
				if f.endswith(nSeq + '.trees.txt') and f.startswith(('initial_', 'sim')):
					logFiles.append(f)
					directories.append(root)

# Get summary tables for each Beast run
def helperF(f,d):
	return(treePosteriorDistributionStatistics(f, d, beastModel, regions))

with ProcessPoolExecutor(max_workers=10) as executor:
	
	lineage = []
	introduction = []

	for result in executor.map(helperF, logFiles, directories):
		lineage.append(result[0])
		introduction.append(result[1])
        
        
# Concatenate dataframes
lineage_pd = pd.concat(lineage, axis = 0).reset_index(drop = True)
introduction_pd = pd.concat(introduction, axis = 0).reset_index(drop = True)

# Write the dataframe
if beastModel == "files": beastModel = 'sim'
lineage_pd.to_csv(directory + 'analyses/lineage_migration_' + beastModel + '_' + nSeq + '.txt', 
	index = False, header = True, sep = '\t', na_rep = 'NA')
introduction_pd.to_csv(directory + 'analyses/introduction_dates_' + beastModel + '_' + nSeq + '.txt', 
	index = False, header = True, sep = '\t', na_rep = 'NA')