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
sys.path.insert(0, os.path.abspath('Python/'))
from treeTopologies import *

# Get arguments
simulationFramework = sys.argv[1]   # "7demes" or "3demes" 
models = sys.argv[2:]               # "dta", "basta", "mascot" or "glm"

# Change 
directory = simulationFramework
os.chdir(directory)

# List of files
treeList = []
dirList = []

for root, subdirs, files in os.walk('.'):
    if root.endswith('dta') :
        for f in files:
            if '.mcc.tree' in f:
                treeList.append(f)
                dirList.append(root)

# Helper Function
def helperF(t, d):
    return(compareBeasttoSimulation(t, simulationFramework, d, *models))

# Perform linear regressions
with ProcessPoolExecutor(max_workers=10) as executor : 
    out = []
    for result in executor.map(helperF, treeList, dirList):
        out.append(result)

# Write results
out = pd.concat(out, sort = False, ignore_index = True)
out.to_csv('analyses/topology_regression.txt', 
	index = False, header = True, sep = "\t", na_rep = 'NA')