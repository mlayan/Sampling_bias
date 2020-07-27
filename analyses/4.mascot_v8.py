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
import multiprocessing as mp
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor

# Absolute paths
#stAbs = "/mnt/gaia/"  # Local computer
stAbs = "/pasteur/projets/common/"  # Work on the cluster

## Import custom functions
sys.path.append(stAbs + 'MMMI_Rage/Python_modules')
from geneticAndSpatialEstimates import *

# runType
runType = "mascot_v8"


##########################################
## Beast 2 - Mascot - M1
##########################################
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

print("\nM1 computational time:")
print(end-start)

# Concatenate dataframes
data = pd.concat(out, axis = 0)

# Write the dataframe
data.to_csv('analyses/log_files_' + runType + '_ESS.txt', 
            index = False, 
            header = True, 
            sep = '\t', 
            na_rep = 'NA')