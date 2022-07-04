#!/usr/bin/env python
# coding: utf-8

"""
CODE TO ANALYZE BEAST CHAINS ON RABV
"""

# Import libraries
from os import listdir
from os.path import isfile, join, abspath
import sys
import re
import pandas as pd
import numpy as np
import dendropy as dp
from pathlib import Path
from lxml import *
from io import StringIO
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import datetime, timedelta
from math import ceil
import matplotlib.pyplot as plt

# Import local modules
sys.path.append('python/')
from mathFunctions import * 
from geneticAndSpatialEstimates import *

## Change directory
os.chdir('rabies/')
regions = ["Luzon", "OrientalMindoro", "Cebu", "NegrosOriental", "Mindanao", "Catanduanes"]
geneticParams = ['gammaShape','kappa','ucldMean','ucldStdev','rate_coefficientOfVariation']
burnIn = 0.1
allresults = list()
mrst=2010.123


#####################################################
# Load DTA
#####################################################
allDta = list()
for d in ["dta1", "dta2", "dta3"]:
    temp_dta = pd.read_csv("dta/tohma_" + d + ".log.txt", sep = "\t", comment = "#")
    b = math.ceil(temp_dta.shape[0] * burnIn)
    temp_dta = temp_dta.iloc[b:temp_dta.shape[0]]
    allDta.append(temp_dta)
dta = pd.concat(allDta).reset_index(drop=True)
stepInterval = dta['state'].iloc[1] - dta['state'].iloc[0]

# Rename columns
dta = renameColumns(dta, "dta", "tohma")
dta.rename(columns={
    'alpha':'gammaShape',
    'default.ucld.mean':'ucldMean',
    'default.ucld.stdev':'ucldStdev',
    'default.covariance':'rate_covariance',
    'default.coefficientOfVariation':'rate_coefficientOfVariation'
}, inplace = True)
dta.columns = dta.columns.str.replace("c_counts", "nMigration").str.replace("\[1\]|c_", '', regex = True)
dta.columns = dta.columns.str.replace("traits.rates.", "forwards_", regex = False)
dta.columns = dta.columns.str.replace("traits.indicators.", "indicators_", regex = False)
dta.columns = dta.columns.str.replace("\.", "_", regex = True)
dta.columns = dta.columns.str.replace(" ", "", regex = True)

# Write total migration events to file
posteriorDistnMig(dta, "dta").to_csv("figures/distribution_dta.txt", sep="\t", index=False)

# Write root age distribution
rootAge = dta[["state", "ageRoot"]].copy()
rootAge.state = "dta"
rootAge.rename(columns={"state":"model"}).to_csv("figures/ageroot_dta.txt", sep ="\t", index=False)

# Write genetic parameter distribution
geneticParams_df = dta[geneticParams + ["rate_covariance"]].copy()
geneticParams_df["model"] = "dta"
pd.melt(geneticParams_df, id_vars=["model"], value_vars=geneticParams + ["rate_covariance"], 
        var_name="parameter", value_name='value').to_csv("figures/geneticparams_dta.txt", sep ="\t", index=False)

# Dataframe with the root location
dtaRootLocation = dta[['traits']].reset_index(drop=True)
dtaRootLocation = dtaRootLocation.rename(columns={"traits":"rootLocation"})
dta.drop(columns=['traits'], inplace = True)

# Summary statistics
all_percentiles = [x for x in range(5, 100, 5)] + [1, 2.5, 99, 97.5]
all_percentiles = sorted([x/100 for x in all_percentiles])
results = dta.drop(columns=['state']).describe(percentiles = all_percentiles).transpose(copy = True)
#results = dta.drop(columns=['state']).describe(percentiles = [0.025, 0.5, 0.975]).transpose(copy = True)
results.drop(['count'], axis = 1, inplace = True)

# ESS of all continuous parameters 
# Migration rates are excluded when BSSVS is implemented
discrete_or_rates = [x for x in dta.columns if x in ["state", "allTransitions", "traits_nonZeroRates", "traits_count"]] + [x for x in dta.columns if "indicators" in x or "nMigration" in x or "forwards" in x]
results['ESS'] = dta.drop(columns=discrete_or_rates).apply(
    lambda x: ESS(x, stepInterval)
    )

# 95% HPD intervals 
results['2_5_hpd'] = dta.drop(columns=['state']).apply(
    lambda x: hpd(x, lower_only = True)
    ) 
results['97_5_hpd'] = dta.drop(columns=['state']).apply(
    lambda x: hpd(x, upper_only = True)
    )

# Change results index
results['parameter'] = results.index
results.reset_index(inplace = True, drop = True)

##############################################
# Deal with spatial data when BSSVS is 
# implemented and compute conditional ess  
##############################################
bssvsResults = bssvsStatistics(dta, len(regions), False, "dta", stepInterval, "true")
results = pd.concat([results, bssvsResults], sort = False, ignore_index = True) 

# Remove lines corresponding to raw migration rates
results = results[~results.parameter.str.match(r'^backwards_.*$|^forwards_.*$|^nMigration_.*$')]
results.reset_index(inplace = True, drop = True)

##############################################
# Add root location if wanted
##############################################
rootResults = rootLocationAnalysis(None, "dta/", "dta", regions, dtaRootLocation, dict(zip(regions, regions)), 10)

results = pd.concat([results, rootResults], sort = False, ignore_index = True)
results.rename(columns={"2.5%":"2_5", "97.5%":"97_5", "50%":"50"}, inplace = True)
results.columns = results.columns.str.replace("%", "")
results["model"] = "dta"
results.to_csv("figures/output_dta.txt", sep="\t", index=False)





#####################################################
# Load BASTA
#####################################################
allBasta = list()

for d in ["basta1", "basta2", "basta3"]:
    temp_basta = pd.read_csv("basta/tohma_" + d + ".log.txt", sep = "\t", comment = "#")
    b = math.ceil(temp_basta.shape[0] * burnIn)
    temp_basta = temp_basta.iloc[b:temp_basta.shape[0]]
    allBasta.append(temp_basta)

basta = pd.concat(allBasta).reset_index(drop=True)
stepInterval = basta['Sample'].iloc[1] - basta['Sample'].iloc[0]

# Remove burn in and rename columns
basta.rename(columns={
    'gammaShape.s:tohma_sequences_aln':'gammaShape',
    'kappa.s:tohma_sequences_aln': "kappa",
    'ucldMean.c:tohma_sequences_aln':'ucldMean',
    'ucldStdev.c:tohma_sequences_aln':'ucldStdev',
    'rate.c:tohma_sequences_aln.mean':'rate_mean',
    'rate.c:tohma_sequences_aln.variance':'rate_variance',
    'rate.c:tohma_sequences_aln.coefficientOfVariation':'rate_coefficientOfVariation'
}, inplace = True)
bastaRegionDict = {"0": "Catanduanes", "1":"Cebu", "2":"Luzon", "3":"Mindanao", "4":"NegrosOriental", "5":"OrientalMindoro"}
basta = renameColumns(basta, "basta", "tohma", bastaRegionDict)
basta = forwardsRates(basta, "forwards", "basta")

# Write total migration events to file
posteriorDistnMig(basta, "basta").to_csv("figures/distribution_basta.txt", 
                                         sep="\t", index=False)

# Write genetic parameter distribution
geneticParams_df = basta[geneticParams + ["rate_variance"]].copy()
geneticParams_df["model"] = "basta"
pd.melt(geneticParams_df, id_vars=["model"], value_vars=geneticParams + ["rate_variance"], 
        var_name="parameter", value_name='value').to_csv("figures/geneticparams_basta.txt",
                                                         sep ="\t", index=False)

# Dataframe with the root location
bastaRootLocation = basta[['rootLocation']].reset_index(drop=True)
basta.drop(columns = ["rootLocation"], inplace = True)

# Summary statistics
results = basta.drop(columns=['state']).describe(percentiles = [0.025, 0.5, 0.975]).transpose(
    copy = True)
results.drop(['count'], axis = 1, inplace = True)

# ESS of all continuous parameters 
# Migration rates are excluded when BSSVS is implemented
discrete_or_rates = [x for x in basta.columns if x in ["state"]] + [x for x in basta.columns if "indicators" in x or 
                     "nMigration" in x or 'backwards' in x]
results['ESS'] = basta.drop(columns=discrete_or_rates).apply(
    lambda x: ESS(x, stepInterval)
    )

# 95% HPD intervals 
results['2_5_hpd'] = basta.drop(columns=['state']).apply(
    lambda x: hpd(x, lower_only = True)
    ) 
results['97_5_hpd'] = basta.drop(columns=['state']).apply(
    lambda x: hpd(x, upper_only = True)
    )

# Change results index
results['parameter'] = results.index
results.reset_index(inplace = True, drop = True)

##############################################
# Deal with spatial data when BSSVS is 
# implemented and compute conditional ess  
##############################################
bssvsResults = bssvsStatistics(basta, len(regions), False, "basta", stepInterval, "all")
results = pd.concat([results, bssvsResults], sort = False, ignore_index = True) 

# Remove lines corresponding to raw migration rates
results = results[~results.parameter.str.match(r'^backwards_.*$|^forwards_.*$|^nMigration_.*$')]
results.reset_index(inplace = True, drop = True)

##############################################
# Add root location if wanted
##############################################
rootResults = rootLocationAnalysis(None, "basta/", "basta", regions, bastaRootLocation, 
                                   dict(zip(regions, regions)), 10)
results = pd.concat([results, rootResults], sort = False, ignore_index = True)
results.rename(columns={"2.5%":"2_5", "97.5%":"97_5", "50%":"50"}, inplace = True)
results["model"]= "basta"
results.to_csv("figures/output_basta.txt", sep="\t", index=False)


#####################################################
# Merge basta trees from posterior distribution
#####################################################
treeList = dp.TreeList()
ageRoot = []
files = [
        "basta/tohma_basta1.trees.txt", 
        "basta/tohma_basta2.trees.txt",
        "basta/tohma_basta3.trees.txt"
    ]

for f in files:
    treeY = dp.Tree.yield_from_files(
        files = [f],
        schema = 'nexus',  
        preserve_underscores = True,  
        extract_comment_metadata = True,
        rooting = 'default-rooted'
    )
    b = 101
    
    for i, tree in enumerate(treeY):
        if i >= b:
            ageRoot.append(mrst - max(tree.calc_node_ages(is_force_max_age =True)))
            treeList.append(tree)

treeList.write(path="basta/tohma_basta_all.trees.txt", schema="nexus")
pd.DataFrame({"model":"basta","ageRoot":ageRoot}).to_csv("figures/ageroot_basta.txt", sep ="\t", index=False)



#####################################################
# Load MASCOT
#####################################################
allMascot = list()
mascotFiles = ["mascot1", "mascot2", "mascot3"]

for d in mascotFiles:
    temp_mascot = pd.read_csv("mascot/tohma_" + d + ".log.txt", sep = "\t", comment = "#")
    b = math.ceil(temp_mascot.shape[0] * burnIn)
    temp_mascot = temp_mascot.iloc[b:temp_mascot.shape[0]]
    allMascot.append(temp_mascot)

mascot = pd.concat(allMascot).reset_index(drop=True)
mascot.rename(columns={"proportionInvariant":"pInv", 
                       "Sample":"state",
                       "treeLikelihood.tohma_sequences_aln": "treeLikelihood",
                       "sum(indicatorsConstantBSSVS.)": "sumNonZeroRates"
                      }, inplace=True)
stepInterval = mascot['state'].iloc[1] - mascot['state'].iloc[0]

# Remove burn in and rename columns
mascot.columns = mascot.columns.str.replace("Negros_Oriental", "NegrosOriental", regex = False)
mascot.columns = mascot.columns.str.replace("Oriental_Mindoro", "OrientalMindoro", regex = False)
mascot.columns = mascot.columns.str.replace("migrationEvents", "nMigration", regex = False)
mascot.columns = mascot.columns.str.replace("b_migration", "backwards", regex = False)
mascot.columns = mascot.columns.str.replace("to_", "", regex = False)
mascot.columns = mascot.columns.str.replace(".", "_", regex = False)
mascot = forwardsRates(mascot, "forwards", "mascot")

# Write root age distribution
rootAge = mascot[["state", "TreeHeight"]].copy()
rootAge.TreeHeight = mrst - rootAge.TreeHeight
rootAge.state = "mascot"
rootAge.rename(columns={"state":"model", "TreeHeight":"ageRoot"}).to_csv("figures/ageroot_mascot.txt", sep ="\t", index=False)

# Write genetic parameter distribution
geneticParams_df = mascot[geneticParams + ["rate_variance"]].copy()
geneticParams_df["model"] = "mascot"
pd.melt(geneticParams_df, id_vars=["model"], value_vars=geneticParams + ["rate_variance"], 
        var_name="parameter", value_name='value').to_csv(
    "figures/geneticparams_mascot.txt", sep ="\t", index=False)

# Write total migration events to file
posteriorDistnMig(mascot, "mascot").to_csv("figures/distribution_mascot.txt", 
                                           sep="\t", index=False)


# Summary statistics
results = mascot.drop(columns=['state']).describe(percentiles = [0.025, 0.5, 0.975]).transpose(copy = True)
results.drop(['count'], axis = 1, inplace = True)

# ESS of all continuous parameters 
# Migration rates are excluded when BSSVS is implemented
discrete_or_rates = [x for x in mascot.columns if x in ["state", "sumNonZeroRates", "nrMigrationEvents"]] + [x for x in mascot.columns if "nMigration" in x or 'backwards' in x]
results['ESS'] = mascot.drop(columns=discrete_or_rates).apply(
    lambda x: ESS(x, stepInterval)
    )

# 95% HPD intervals 
results['2_5_hpd'] = mascot.drop(columns=['state']).apply(
    lambda x: hpd(x, lower_only = True)
    ) 
results['97_5_hpd'] = mascot.drop(columns=['state']).apply(
    lambda x: hpd(x, upper_only = True)
    )

# Change results index
results['parameter'] = results.index
results.reset_index(inplace = True, drop = True)

##############################################
# Deal with spatial data when BSSVS is 
# implemented and compute conditional ess  
##############################################
bssvsResults = bssvsStatistics(mascot, len(regions), False, "mascot", stepInterval, "all")
results = pd.concat([results, bssvsResults], sort = False, ignore_index = True) 

# Remove lines corresponding to raw migration rates
results = results[~results.parameter.str.match(r'^backwards_.*$|^forwards_.*$|^nMigration_.*$')]
results.reset_index(inplace = True, drop = True)

# Root location
ann = annotationDict["mascot"]
rootLocation = []
allTrees = dp.TreeList()

for x in mascotFiles:
    treeCollectionYielder = dp.Tree.yield_from_files(
        files = ["mascot/tohma_" + x + ".trees.txt"], 
        schema = 'nexus',  
        preserve_underscores = True,  
        extract_comment_metadata = True
    )
    
    for i, tree in enumerate(treeCollectionYielder):
        if i >= b:
            allTrees.append(removeUnifurcations(tree))
            rootLocation.append(tree.seed_node.annotations.get_value(name=ann))

allTrees.write(path="mascot/tohma_mascot_full.trees.txt", schema="nexus")
resultsRoot = pd.DataFrame(pd.Series(rootLocation).value_counts() / len(rootLocation) )
resultsRoot.reset_index(inplace=True)
resultsRoot.columns = ["parameter", "value"]

if resultsRoot.shape[0] != len(regions):
    for x in regions:
        if x not in list(resultsRoot.parameter):
            resultsRoot = resultsRoot.append(pd.DataFrame({'parameter':[x], 'value':[0.0]}), ignore_index = True)

resultsRoot.parameter = 'root_' + resultsRoot.parameter.astype(str)
kl = KLRootPrediction(resultsRoot.value, len(regions), True)
resultsRoot = resultsRoot.append(pd.DataFrame({"parameter": ["KL"], "value": [kl]}), ignore_index = True)

results = pd.concat([results, resultsRoot], sort = False, ignore_index = True)
results.rename(columns={"2.5%":"2_5", "97.5%":"97_5", "50%":"50"}, inplace = True)
results["model"]= "mascot"
results.to_csv("figures/output_mascot.txt", sep="\t", index=False)

