#!/usr/bin/env python
# coding: utf-8

"""
ANALYZE BEAST LOG FILES
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-02-28' 
__last_update__ = '2020-05-13'

# Import libraries
import os
import sys
import re
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

## Import custom functions
sys.path.append('/python/')
from geneticAndSpatialEstimates import *

###############################################
## SUM UP  
# Runtype
runType = "dta"

## Directory
cond = 'mig' + str(sys.argv[1])
directory = "../" + cond
os.chdir(directory)

# Dictionary of the regions
regionDict = {
    1: "Region1",
    2: "Region2",
    3: "Region3",
    4: "Region4",
    5: "Region5",
    6: "Region6",
    7: "Region7"
}

# Load log files
logFiles = []
directories = []

for root, subdirs, files in os.walk('.'):
	if runType in root:
		for f in files:
			if 'log.txt' in f:
				logFiles.append(f)
				directories.append(root)


# Get summary tables for each Beast run
def helperF(f,d):
    return(logFileWrangler(f, d, cond, runType, regionDict,
        rootLocation = True, ess = True))

with ProcessPoolExecutor() as executor:
    out = []
    
    for result in executor.map(helperF, logFiles, directories):
        out.append(result)

# Concatenate dataframes
data = pd.concat(out, axis = 0)

# Write the dataframe
data.to_csv('../analyses/log_files_' + runType + '_' + cond + '_ESS.txt', 
            index = False, 
            header = True, 
            sep = '\t', 
            na_rep = 'NA')
