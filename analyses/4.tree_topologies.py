#!/usr/bin/env python
# coding: utf-8

"""
COMPUTE THE REGRESSION COEFFICIENT COEFFICIENT 
OF THE PAIRWISE tMRCA
"""

__author__ = "Maylis Layan"
__creation_date__ = "2020/05/25"
__last_update__ = "2020/05/26"


# Import modules
import os
import sys
import re
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

# Import personnal modules 
sys.path.insert(0, os.path.abspath('python'))
from treeTopologies import *


# Change working directory
directory = absPath + "HKY_radiation"
matrix = "radiation"
#int(re.match(r'.*M(\d+)', directory).group(1))
models = 'dta'

# List of files
treeList = []
dirList = []

for entry in os.scandir():
    if entry.is_dir() and 'simulation' in entry.name:
        
        nSim = int(re.match(r'simulation(\d+)$', entry.name).group(1))
        if nSim < 11:
            
            for files in os.scandir(entry.name + '/dta') : 
                if files.is_file() and '.mcc.tree' in files.name:
                    treeList.append(files.name)
                    dirList.append(entry.name + '/dta')

# Helper Function
def helperF(t, d):
    return(compareBeasttoSimulation(t, matrix, d, *models))

# Perform linear regressions
with ProcessPoolExecutor() as executor : 
    out = []
    for result in executor.map(helperF, treeList, dirList):
        out.append(result)

# Write results
out = pd.concat(out, sort = False, ignore_index = True)
out.to_csv('inputfiles/topology_regression_sim1to10.txt', index = False, 
    header = True, sep = "\t", na_rep = 'NA')

