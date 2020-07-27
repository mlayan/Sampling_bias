#!/usr/bin/env python
# coding: utf-8

#!/usr/bin/python
# -*- coding: utf-8 -*-

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
import math
import time
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

## Import custom functions
sys.path.append('python')
from geneticAndSpatialEstimates import *

# runType
runType = "mascot_v8"

## Directory
directory = stAbs + "MMMI_Rage/HKY_M1/"
os.chdir(directory)

# Matrix
nMatrix = int(re.search(r'M([\d]+)/', directory).group(1))

# Load log files
logFiles = []
directories = []

for entry in os.scandir():    
    if entry.is_dir() and entry.name.startswith('simulation'):
        nSim = int(re.search(r'\d+$', entry.name).group(0))
        
        if nSim < 11:
            for files in os.scandir(entry.name + "/" + runType):
                if files.is_file() and files.name.endswith('log.txt') \
                						and 'stratified' not in files.name :
                    logFiles.append(files.name)
                    directories.append(entry.name + "/" + runType)


# Get summary tables for each Beast run
def helperF(f,d):
    return(logFileWrangler(f, d, nMatrix, runType, 
        rootLocation = True, forwards = "all", ess = True))

start = time.time()

with ProcessPoolExecutor() as executor:
    out = []
    
    for result in executor.map(helperF, logFiles, directories):
        out.append(result)
        
end = time.time()

print("\nComputational time:")
print(end-start)
print("\n")

# Concatenate dataframes
data = pd.concat(out, axis = 0)

# Write the dataframe
data.to_csv('inputfiles/log_files_' + runType + '_ESS.txt', 
            index = False, 
            header = True, 
            sep = '\t', 
            na_rep = 'NA')
