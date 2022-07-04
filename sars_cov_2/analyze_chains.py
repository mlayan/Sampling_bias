#!/usr/bin/env python
# coding: utf-8

"""
CODE TO ANALYZE BEAST CHAINS ON SARS-COV-2 
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
os.chdir('sars_cov_2/')
regions = ["Africa", "Americas", "Asia", "China", "Europe", "Oceania"]
geneticParams = ['kappa','gammaShape','pInv','clockRate']
burnIn = 0.1
mrst = 2020.175
allresults = list()


#####################################################
# Load DTA
#####################################################
allDta = list()
for d in ["dta", "dta2", "dta3"]:
    temp_dta = pd.read_csv("dta/sars_cov_2_" + d + ".log.txt", sep = "\t", comment = "#")
    b = math.ceil(temp_dta.shape[0] * burnIn)
    temp_dta = temp_dta.iloc[b:temp_dta.shape[0]]
    allDta.append(temp_dta)
dta = pd.concat(allDta).reset_index(drop=True)
stepInterval = dta['state'].iloc[1] - dta['state'].iloc[0]

# Rename columns
dta = renameColumns(dta, "dta", "sars_cov_2")
dta.rename(columns = {"alpha":"gammaShape"}, inplace=True)
dta.columns = dta.columns.str.replace("c_counts", "nMigration").str.replace("\[1\]|c_", '', regex = True)
dta.columns = dta.columns.str.replace("continent.rates.", "forwards_", regex = False)
dta.columns = dta.columns.str.replace("continent.indicators.", "indicators_", regex = False)
dta.columns = dta.columns.str.replace("\.", "_", regex = True)

# Write root age distribution
rootAge = dta[["state", "ageRoot"]].copy()
rootAge.state = "dta"
rootAge.rename(columns={"state":"model"}).to_csv("figures/ageroot_dta.txt", sep ="\t", index=False)

# Write genetic parameter distribution
geneticParams_df = dta[geneticParams].copy()
geneticParams_df["model"] = "dta"
pd.melt(geneticParams_df, id_vars=["model"], value_vars=geneticParams, 
        var_name="parameter", value_name='value').to_csv("figures/geneticparams_dta.txt", sep ="\t", index=False)


# Write total migration events to file
posteriorDistnMig(dta, "dta").to_csv("figures/distribution_dta.txt", sep="\t", index=False)

# Dataframe with the root location
dtaRootLocation = dta[['continent']].reset_index(drop=True)
dtaRootLocation = dtaRootLocation.rename(columns={"continent":"rootLocation"})
dta.drop(columns=['continent'], inplace = True)

# Summary statistics
results = dta.drop(columns=['state']).describe(percentiles = [0.025, 0.5, 0.975]).transpose(copy = True)
results.drop(['count'], axis = 1, inplace = True)

# ESS of all continuous parameters 
# Migration rates are excluded when BSSVS is implemented
discrete_or_rates = [x for x in dta.columns if x in ["state", "allTransitions", "continent_nonZeroRates", "continent_count"]] + [x for x in dta.columns if "indicators" in x or "nMigration" in x or "forwards" in x]
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

##############################################
# verify ESS values
##############################################
nLow = sum(~results.ESS.isna() & results.ESS < 200 & (results.BF >=3 | results.BF.isna()))
if nLow:
    print(results[~results.ESS.isna() & results.ESS < 200 & (results.BF >=3 | results.BF.isna())])

results["model"] = "dta"
results.to_csv("figures/output_dta.txt", sep="\t", index=False)
dta.to_csv("figures/posterior_dta.txt", sep ="\t", index=False)







#####################################################
# Load BASTA
#####################################################
allBasta = list()

for d in ["basta1", "basta2", "basta3"]:
    temp_basta = pd.read_csv("basta/sars_cov_2_" + d + ".log.txt", sep = "\t", comment = "#")
    b = math.ceil(temp_basta.shape[0] * burnIn)
    temp_basta = temp_basta.iloc[b:temp_basta.shape[0]]
    allBasta.append(temp_basta)

basta = pd.concat(allBasta).reset_index(drop=True)
stepInterval = basta['Sample'].iloc[1] - basta['Sample'].iloc[0]

# Remove burn in and rename columns
basta.columns = basta.columns.str.replace("proportionInvariant.s:sars_cov_2", "pInv", regex=False)
basta.columns = basta.columns.str.replace("gammaShape.s:sars_cov_2", "gammaShape", regex = False)
bastaRegionDict = {"0":"Africa", "1":"Americas", "2":"Asia", "3": "China", "4":"Europe", "5":"Oceania"}
basta = renameColumns(basta, "basta", "sars_cov_2", bastaRegionDict)
basta = forwardsRates(basta, "forwards", "basta")

# Write genetic parameter distribution
geneticParams_df = basta[geneticParams].copy()
geneticParams_df["model"] = "basta"
pd.melt(geneticParams_df, id_vars=["model"], value_vars=geneticParams, 
        var_name="parameter", value_name='value').to_csv("figures/geneticparams_basta.txt", sep ="\t", index=False)
# Write total migration events to file
posteriorDistnMig(basta, "basta").to_csv("figures/distribution_basta.txt", sep="\t", index=False)

# Dataframe with the root location
bastaRootLocation = basta[['rootLocation']].reset_index(drop=True)
basta.drop(columns = ["rootLocation"], inplace = True)

# Sumamry statistics
results = basta.drop(columns=['state']).describe(percentiles = [0.025, 0.5, 0.975]).transpose(copy = True)
results.drop(['count'], axis = 1, inplace = True)

# ESS of all continuous parameters 
# Migration rates are excluded when BSSVS is implemented
discrete_or_rates = [x for x in basta.columns if x in ["state"]] + [x for x in basta.columns if "indicators" in x or "nMigration" in x or 'backwards' in x]
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
rootResults = rootLocationAnalysis(None, "basta/", "basta", regions, bastaRootLocation, dict(zip(regions, regions)), 10)
results = pd.concat([results, rootResults], sort = False, ignore_index = True)
results.rename(columns={"2.5%":"2_5", "97.5%":"97_5", "50%":"50"}, inplace = True)

##############################################
# verify ESS values
##############################################
nLow = sum(~results.ESS.isna() & results.ESS < 200 & (results.BF >=3 | results.BF.isna()))
if nLow:
    print(results[~results.ESS.isna() & results.ESS < 200 & (results.BF >=3 | results.BF.isna())])

results["model"]= "basta"
results.to_csv("figures/output_basta.txt", sep="\t", index=False)
basta.to_csv("figures/posterior_basta.txt", sep ="\t", index=False)


#####################################################
# Separate post burn-in trees according to 
# the mode of the tree prior
#####################################################
majMode = dp.TreeList()
minMode = dp.TreeList()

for d in ["basta1", "basta2", "basta3"]:
    temp_basta = pd.read_csv("basta/sars_cov_2_" + d + ".log.txt", sep = "\t", comment = "#")
    b = math.ceil(temp_basta.shape[0] * burnIn)
    
    treeY = dp.Tree.yield_from_files(
        files = ["basta/sars_cov_2_" + d+ ".trees.txt"],
        schema = 'nexus',  
        preserve_underscores = True,  
        extract_comment_metadata = True,
        rooting = 'default-rooted'
    )
    
    for i, tree in enumerate(treeY):
        if i >= b and i < temp_basta.shape[0] and temp_basta["mutationRate"].iloc[i] >= 11e-4:
            majMode.append(tree)
        
        if i >= b and i < temp_basta.shape[0]  and temp_basta["mutationRate"].iloc[i] <= 10e-4:
            minMode.append(tree)

majMode.write(path="basta/sars_cov_2_basta_majmode.trees.txt", schema="nexus")
minMode.write(path="basta/sars_cov_2_basta_minmode.trees.txt", schema="nexus")






#####################################################
# Load MASCOT
#####################################################
allMascot = list()
mascotFiles = ["mascot1", "mascot2", "mascot3"]

for d in mascotFiles:
    temp_mascot = pd.read_csv("mascot/sars_cov_2_" + d + ".log.txt", sep = "\t", comment = "#")
    b = math.ceil(temp_mascot.shape[0] * burnIn)
    temp_mascot = temp_mascot.iloc[b:temp_mascot.shape[0]]
    allMascot.append(temp_mascot)

mascot = pd.concat(allMascot).reset_index(drop=True)
mascot.rename(columns={"proportionInvariant":"pInv", 
                       "Sample":"state",
                       "treeLikelihood.sars_cov_2": "treeLikelihood",
                       "sum(indicatorsConstantBSSVS_)": "sumNonZeroRates"
                      }, inplace=True)
stepInterval = mascot['state'].iloc[1] - mascot['state'].iloc[0]

# Remove burn in and rename columns
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
geneticParams_df = mascot[geneticParams].copy()
geneticParams_df["model"] = "mascot"
pd.melt(geneticParams_df, id_vars=["model"], value_vars=geneticParams, 
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
for x in mascotFiles:
    treeCollectionYielder = dp.Tree.yield_from_files(
        files = ["mascot/sars_cov_2_" + x + ".trees.txt"], 
        schema = 'nexus',  
        preserve_underscores = True,  
        extract_comment_metadata = True
    )
    
    for i, tree in enumerate(treeCollectionYielder):
        if i >= b:
            rootLocation.append(tree.seed_node.annotations.get_value(name=ann))

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

results = pd.concat([results, rootResults], sort = False, ignore_index = True)
results.rename(columns={"2.5%":"2_5", "97.5%":"97_5", "50%":"50"}, inplace = True)

##############################################
# verify ESS values
##############################################
nLow = sum(~results.ESS.isna() & results.ESS < 200 & (results.BF >=3 | results.BF.isna()))
if nLow:
    print(results[~results.ESS.isna() & results.ESS < 200 & (results.BF >=3 | results.BF.isna())])

results["model"]= "mascot"
results.to_csv("figures/output_mascot.txt", sep="\t", index=False)
mascot.to_csv("figures/posterior_mascot.txt", sep ="\t", index=False)


##############################################
# Create full tree posterior distribution
##############################################
treeList = dp.TreeList()

for x in mascotFiles:
    treeList2 = dp.TreeList()
    treeCollectionYielder = dp.Tree.yield_from_files(
        files = ["mascot/sars_cov_2_" + x + ".trees.txt"], 
        schema = 'nexus',  
        preserve_underscores = True,  
        extract_comment_metadata = True,
        rooting = 'default-rooted'
    )
    
    logFile = pd.read_csv("mascot/sars_cov_2_" + x + ".log.txt", comment='#', sep="\t")
    b = math.ceil(logFile.shape[0] * burnIn / 100)
    
    for i, tree in enumerate(treeCollectionYielder):
        if i > b:
            treeList.append(removeUnifurcations(tree))
            
treeList.write(path="mascot/sars_cov_2_mascot_full.trees.txt", schema="nexus")









#####################################################
# Load MASCOT-GLM
#####################################################
# Log file
allresults = []
for g in ["glm_wid"]: 
    allGlm = list()
    glmFiles = [g.replace("glm_", "") + str(i) for i in range(1,4)]

    for d in glmFiles:
        temp_glm = pd.read_csv(g + "/sars_cov_2_" + d + ".log.txt", sep = "\t", comment = "#")
        b = math.ceil(temp_glm.shape[0] * burnIn)
        temp_glm = temp_glm.iloc[b:temp_glm.shape[0]]
        allGlm.append(temp_glm)
    log = pd.concat(allGlm).reset_index(drop=True)

    stepInterval = log['Sample'].iloc[1] - log['Sample'].iloc[0]

    # Remove Ne.Region columns
    toDrop = [x for x in log.columns if "Ne." in x]
    log.drop(columns=toDrop, inplace=True)

    # Rename columns
    log.rename(columns={
        "Sample":"state",
        "sum(NeIndicatorGLM.)":"sumNePredNonZero", 
        "sum(migrationIndicatorGLM.)":"sumMigPredNonZero", 
        "proportionInvariant": "pInv", 
        "treeLikelihood.sars_cov_2": "treeLikelihood"
        }, inplace = True)
    log.columns = log.columns.str.replace("migrationEvents", "nMigration", regex = False)
    log.columns = log.columns.str.replace("mig.", "forwards_", regex = False)
    log.columns = log.columns.str.replace("to_", "", regex = False)
    log.columns = log.columns.str.replace(".", "_", regex = False)
    scalers = log.filter(regex=("scaler_.*"))
    log.drop(columns=[x for x in log.columns if "scaler" in x], inplace = True)
    
    # Write root age distribution
    rootAge = log[["state", "TreeHeight"]].copy()
    rootAge.TreeHeight = mrst - rootAge.TreeHeight
    rootAge.state = g
    rootAge.rename(columns={"state":"model", "TreeHeight":"ageRoot"}).to_csv(
        "figures/ageroot_" + g + ".txt", sep ="\t", index=False)

    # Write genetic parameter distribution
    geneticParams_df = log[geneticParams].copy()
    geneticParams_df["model"] = g
    pd.melt(geneticParams_df, id_vars=["model"], value_vars=geneticParams, 
            var_name="parameter", value_name='value').to_csv(
        "figures/geneticparams_" + g + ".txt", sep ="\t", index=False)

    # Write total migration events to file
    posteriorDistnMig(log, g).to_csv("figures/distribution_" + g+ ".txt", sep="\t", index=False)
    
    ##############################################
    # Summary statistics and ESS of all parameters
    ##############################################
    # Summary statistics (min, mean, median, max, 95% CI)
    results = log.drop(columns=['state']).describe(percentiles = 
        [0.025, 0.5, 0.975]).transpose(copy = True)
    results.drop(['count'], axis = 1, inplace = True)

    # ESS of all continuous parameters 
    # Migration rates are excluded when BSSVS is implemented
    discrete_or_rates = [x for x in log.columns if x in ["state", "sumMigPredNonZero", "sumNePredNonZero", "nrMigrationEvents"]] +                         [x for x in log.columns if "nMigration" in x]

    results['ESS'] = log.drop(columns=discrete_or_rates).apply(
        lambda x: ESS(x, stepInterval)
        )

    # 95% HPD intervals 
    results['2_5_hpd'] = log.drop(columns=['state']).apply(
        lambda x: hpd(x, lower_only = True)
        ) 
    results['97_5_hpd'] = log.drop(columns=['state']).apply(
        lambda x: hpd(x, upper_only = True)
        )

    # Change results index
    results['parameter'] = results.index
    results.reset_index(inplace = True, drop = True)

    # Inclusion of probabilities derived from scalers
    results_scalers = scalers.apply(lambda x: 1-mean(x==0)).to_frame().reset_index()
    results_scalers.columns = ["parameter", "inclusion_prob"]
    results_scalers.parameter = results_scalers.parameter.str.replace("scaler_", "").replace("Ne", "NeGLM_Clock")

    # Merge sum stats and inclusion probabilities
    results = results.merge(results_scalers, how="outer", on = "parameter")

    ##############################################
    # Add root location and write full 
    # posterior tree distribution without burn-in
    ##############################################
    ann = annotationDict["mascot"]
    rootLocation = []
    treePostDist = dp.TreeList()

    for x in glmFiles:
        treeCollectionYielder = dp.Tree.yield_from_files(
            files = [g + "/sars_cov_2_" + x + ".trees.txt"], 
            schema = 'nexus',  
            preserve_underscores = True,  
            extract_comment_metadata = True
        )

        for i, tree in enumerate(treeCollectionYielder):
            if i >= b:
                rootLocation.append(tree.seed_node.annotations.get_value(name=ann))
                treePostDist.append(removeUnifurcations(tree))

    treePostDist.write(path= g + "/sars_cov_2_" + g + "_full.trees.txt", schema="nexus")
    resultsRoot = pd.DataFrame(pd.Series(rootLocation).value_counts() / len(rootLocation) )
    resultsRoot.reset_index(inplace=True)
    resultsRoot.columns = ["parameter", "value"]

    if resultsRoot.shape[0] != len(regions):
        for x in regions:
            if x not in list(resultsRoot.parameter):
                resultsRoot = resultsRoot.append(pd.DataFrame({'parameter':[x], 'value':[0.0]}), 
                                                 ignore_index = True)

    resultsRoot.parameter = 'root_' + resultsRoot.parameter.astype(str)
    kl = KLRootPrediction(resultsRoot.value, len(regions), True)
    resultsRoot = resultsRoot.append(pd.DataFrame({"parameter": ["KL"], "value": [kl]}), 
                                     ignore_index = True)

    results = pd.concat([results, resultsRoot], sort = False, ignore_index = True)

    results.rename(columns={"2.5%":"2_5", "97.5%":"97_5", "50%":"50"}, inplace = True)
    results["model"]= g
    results.to_csv("figures/output_" + g + ".txt", sep = "\t", index=False)
    log.to_csv("figures/posterior_"+ g +".txt", sep ="\t", index=False)
