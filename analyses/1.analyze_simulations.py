#!/usr/bin/env python
# coding: utf-8

"""
ANALYZE SIMULATION FILES 
- Write Newick trees
- Compute migration events and regions frequencies
- Compute association index of the tree 
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-05-06' 
__last_update__ = '2020-05-'

# Import libraries
import os
import sys
import re
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

## Import custom functions
sys.path.append('/python/')
from simulatedTrees import *

## Directory
cond = 'mig' + str(sys.argv[1])
directory = "../" + cond
os.chdir(directory)

## Region dictionary
regionDic = {
1: "Region1",
2: "Region2",
3: "Region3",
4: "Region4",
5: "Region5",
6: "Region6",
7: "Region7"
}

##########################################
# Load log files
logFiles = []
directories = []
nFiles = 0

for root, subdirs, files in os.walk('.'):
	if 'files' in root:
		for f in files:
			if 'transmission_chain' in f:
				logFiles += [f] * len(protocols)
				directories += [root] * len(protocols)
				nFiles += 1

protocols += protocols * (nFiles - 1)


# Helper function to pass correctly arguments to the mapper
def helperF(f,p,d):
    return(migrationEvents(f, cond, p, regionDic = regionDic,
    	extractNewickTree=True, directory=d))


# Get summary tables for each Beast run
with ProcessPoolExecutor() as executor:
    results = []
    for result in executor.map(helperF, logFiles, protocols, directories):
        results.append(result)

# Concatenate dataframes
data = pd.concat(results, axis = 0)

# Write the dataframe
data.to_csv('../analyses/migration_events_sim1to10_' + cond + '.txt', 
            index = False, 
            header = True, 
            sep = '\t', 
            na_rep = 'NA')
