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
import time
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

## Import custom functions
sys.path.append('python')
from glmEstimates import glmMarkovJumps

# Retrieve arguments 
nDemes = sys.argv[1]            # "3demes" or "7demes"
nSeq = sys.argv[2]              # "150" or "500"
startSim = int(sys.argv[3])     # "1", "11", "21", "31", or "41"
endSim = int(sys.argv[4])       # "10", "20", "30", "40", or "50"
beastModel = "glm"              

# Change working directory
directory = nDemes + "/"
os.chdir(directory)

# Region dictionary
if "3demes" in nDemes:
    regionDict = {
    0:'Region1',
	1:'Region2',
	2:'Region3'
    }
else:
    regionDict = {
    0:'Region1',
    1:'Region2',
    2:'Region3', 
    3:'Region4',
    4:'Region5',
    5:'Region6',
    6:'Region7'
    }

##########################################
# Load log files
logFiles = []
directories = []
simulationsToAnalyze = [ 'simulation' + str(k) for k in range(startSim, endSim + 1) ]


for entry in os.scandir():    
    checkSimulation = [entry.name.endswith(x) for x in simulationsToAnalyze]
    if entry.is_dir() and any(checkSimulation):
        for files in os.scandir(entry.name + "/" + beastModel):
            if files.is_file() and files.name.endswith('log.txt') \
            					and '_' + nSeq + "." in files.name \
            					and files.name.startswith('sim'):

                logFiles.append(files.name)
                directories.append(entry.name + "/" + beastModel)


# Get summary tables for each Beast run
def helperF(f,d):
    return(glmMarkovJumps(
        f, 
        d, 
        nDemes, 
    	regionDict
        ))

with ProcessPoolExecutor(max_workers=12) as executor:
    out = []
    
    for result in executor.map(helperF, logFiles, directories):
        out.append(result)

# Concatenate dataframes
data = pd.concat(out, axis = 0)

# Write the dataframe
data.to_csv('analyses/glm_files_' + nSeq + '_' + str(startSim) + '.txt', 
            index = False, 
            header = True, 
            sep = '\t', 
            na_rep = 'NA')