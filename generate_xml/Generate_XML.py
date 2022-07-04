#!/usr/bin/env python
# coding: utf-8

"""
CODE TO GENERATE XML FILES FOR A SET OF SIMULATED DATA
"""

# Import libraries
import os
import sys
import re
import pandas as pd
import numpy as np
from pathlib import Path

# Import local modules
sys.path.insert(0, os.path.abspath('python/'))
from xmlGenerator import *


## Get all directories containing simulated data
nDemes = ["7demes", "3demes"]
models = ["glm", "mascot", "dta", "basta"]
toAnalyze = [
    "uniformS_10",
    "maxPerRegion_10",  
    "maxPerRegionYear_10", 
    "uniformS_20", 
    "maxPerRegion_20",  
    "maxPerRegionYear_20", 
    "biased_2.5",
    "biased_5",
    "biased_10",
    "biased_20",
    "biased_50",
    "uniform_"
]

listFiles = []
listDirInput = [] 
listDirOutput = []

for md in models: 
    for nd in nDemes:
        for root, subdirs, files in os.walk(nd): 
            if 'files' in root :
                for f in files: 
                    if any([x in f for x in toAnalyze]) and f.endswith(".fasta"):
                        # Location of fasta and trait files
                        dirInput = re.sub(r'^\./', '', root)

                        # Location of XML files to generate
                        dirOutput = dirInput.replace('files', md)
                        testDirOuput = root.replace('files', md)

                        # Lists of directories and file names
                        listDirInput.append(dirInput)
                        listDirOutput.append(dirOutput)
                        listFiles.append(f)


##################################################
# For each simulation produce the XML files 
# for all the sampling strategies and the  
# associated sh files to run BEAST
# Launch them in the meantime
##################################################
for i, f in enumerate(listFiles):

    equalDemes = False
    markovJumps = False
    adjustedBF = None

    nSim = re.match("sim([0-9]+)_", f)
    nSim = int(nSim.group(1))
    
    # Get trait file name
    regionsFile = re.sub('sequences', 'traits', f)
    regionsFile = re.sub('fasta', 'txt', regionsFile)   

    # Logging frequency and chain length
    if "150" in f:
        logEvery = 20000
        chainLength = 20000000
    else:
        logEvery = 40000
        chainLength = 40000000

    # XML template
    if listDirInput[i].lower().endswith("dta"):
        seedValue = None 
        templateFile = "template_beast1_hky_vt.xml"


    if listDirInput[i].lower().endswith("mascot"):
        seedValue = None 
        equalDemes = True
        templateFile = "template_mascot_mapper_vf.xml"
        
    if listDirInput[i].lower().endswith('basta'):
        seedValue = None
        equalDemes = True
        templateFile = "template_basta_vt2.xml"

        
    if listDirOutput[i].lower().endswith("glm"):
        templateFile = "template_mascot_glm_v0.xml" 
        markovJumps = True
        equalDemes = False
        adjustedBF = None
        seedValue = None 

    # Generate the xml file
    create_beast_xml(
        f, 
        regionsFile, 
        templateFile, 
        logEvery, 
        chainLength, 
        dirInput = listDirInput[i], 
        dirOutput = listDirOutput[i], 
        outputName = None, 
        prefix = None, 
        fixedTreeFile = None, 
        fixedTreeDirectory = None, 
        BSSVS = True, 
        markovJumps = markovJumps, 
        equalDemes = equalDemes, 
        adjustedBF = adjustedBF, 
        seed = seedValue
    )


print("Done")

