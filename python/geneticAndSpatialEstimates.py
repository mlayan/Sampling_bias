#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MODULE TO ANALYZE BEAST LOG FILES
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-02-28' 
__last_update__ = '2020-05-04'

# Import libraries
import os
import sys
import re
import math
import pandas as pd
import numpy as np
import dendropy as dp
from itertools import takewhile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
#from scipy.stats import kstest
#from arviz import ess

## Import custom functions
from simulatedTrees import getAncestors
from handleDirectories import checkDirectory
from mathFunctions import *
from mascotMccTrees import removeUnifurcations

## Import regions and log files columns dictionaries
from beastLogDictionaries import * 

# Do not iterate over nodes with a single child (applicable to mascot)
# except if it is the root 
filter_fn = lambda n: len(n.child_nodes()) != 1 or n.parent_node == None  





def logFileWrangler(
	fileName, 
	directory, 
	nMatrix, 
	runType, 
	regionDict,
	bssvs = False,
	rootLocation = False, 
	forwards = "false", 
	burnIn = 10
	):
	"""
	Function to retrieve genetic and spatial estimates from Beast Log Files
		- fileName (dataframe): name of the beast log file  
		- directory (str): directory from which to read input files
		- nMatrix (int): number of the spatial coupling matrix
		- runType (str): type of beast run (beast1, dta, mascot, basta)
		- regionDict (dict): dictionary whose keys are integers and values are the corresponding 
		location strings 
		- rootLocation (bool): if True, the rootLocation probabilities are analyzed
		- forwards (str): 
			argument for mascot beast files only
			if "true", then backwards-in-time migration rates are converted into forwards-in-time migration time 
			if "false", then migration rates are not converted 
			if "all", then both backwards-in-time and forwards-in-time migration rates are analyzed
		- ess (bool): if True, compute ESS for all parameters (depend on the type of beast run) 
			in the Beast log file 

	It returns a dataframe with all genetic and spatial parameters of interest as lines 
	and descritive statistics as columns
	Descriptive statistics:
		- mean
		- std
		- min
		- 2.5-quantile (2.5%) of the credibility interval
		- median (50%)
		- 97.5-quantile (97.5%) of the credibility interval
		- max
		- 2.5-quantile (2.5%_hpd) of the highest posterior density
		- 97.5-quantile (97.5%_hpd) of the highest posterior density
		- BF : Bayes Factor of migration rates when BSSVS is implemented
		- value: only for root location probabilities and Kullback_Leibler divergence of root locations
	Parameters:
		- parameter: name of the parameter
		- nSim: number of the simulation
		- nSeq: sample size
		- protocol: sampling protocol
		- matrix: number of the spatial coupling matrix 
	"""
	print(fileName)

	##############################################
	# Load and rename the dataframe
	##############################################
	# Directory
	directory = checkDirectory(directory)

	# Log file
	if "basta" in runType.lower():
		logs = []
		burnIns = []
		chainLength = []
		for f in [fileName, fileName.replace(".log.txt", "_b.log.txt"), fileName.replace(".log.txt", "_t.log.txt")]:
			if os.path.exists(directory+f):
				log_temp = pd.read_csv(directory+f, sep = "\t", comment = "#")
				chainLength.append(log_temp.shape[0])
				b = round(log_temp.shape[0] * burnIn / 100)
				burnIns.append(b)
				logs.append(log_temp.iloc[b:])
		log = pd.concat(logs, sort = False, ignore_index = True)

	else : 
		logs = []
		for f in [fileName, fileName.replace(".log.txt", "_b.log.txt")]:
			if os.path.exists(directory+f):
				log_temp = pd.read_csv(directory+f, sep = "\t", comment = "#")
				if 'freqParameter.4' in log_temp.columns.tolist() and "mascot" in runType.lower():
					log_temp.rename(columns = {
						'freqParameter.1':'freqParameter.0',
						'freqParameter.2':'freqParameter.1',
						'freqParameter.3':'freqParameter.2',
						'freqParameter.4':'freqParameter.3'
					}, inplace = True)
				logs.append(log_temp)

		log = pd.concat(logs, sort = False, ignore_index = True)
		chainLength = log.shape[0]
		burnIn = round(log.shape[0] * burnIn / 100)
		log = log.iloc[burnIn:]

		
	# Remove spaces and rename columns
	oneRate = [x for x in log.columns if "_and_Region" in x] 
	if len(oneRate):
		j = re.match(r'^b_migration.([^_]*)_.*$', oneRate[0]).group(1)
		i = re.match(r'^.*_([^_]*)$', oneRate[0]).group(1)
		log["b_migration." + i + "_and_" + j ] = log[oneRate[0]]
		log.columns = log.columns.str.replace('_and_Region', '_to_Region')

	# Get region dictionary for Basta chains
	bastaRegionDict = None
	if 'basta' in runType.lower():
		bastaRegionDict = getBastaRegionDict(fileName, directory)

	# Rename columns
	log = renameColumns(log, runType, fileName, bastaRegionDict)
	
	##############################################
	# General informations
	##############################################
	# Information on the data structure
	nSim = int(re.match(r'.*sim(\d*)_', fileName).group(1))
	nSeq = int(re.match(r'.*_(\d*)\.log.*', fileName).group(1))
	protocol = re.match(r'.*sim\d*_(.*)_\d*\.log.*', fileName).group(1)
	
	# Number of regions 
	regions = [re.sub('^.*_', '', x) for x in log.columns.tolist() if 'forwards_' in x or 'backwards_' in x]
	nRegions = len(set(regions))

	# Symmetric or asymmetric matrix 
	if len(regions) == (nRegions * (nRegions-1)) / 2:
		symmetric = True
	elif len(regions) == nRegions * (nRegions - 1):
		symmetric = False
	elif len(regions) == 1 and nRegions == 1:
		symmetric = True
		regions = [x for x in log.columns.tolist() if 'forwards_' in x or 'backwards_' in x][0]
		regions = regions.split("_")
		regions.remove("forwards")
		regions.remove("backwards")
		nRegions = len(regions)
	else:
		raise ValueError("The number of migration rates ({0}) doesn't correspond to \
							the number of regions ({1})".format(len(regions), nRegions))

	# List of regions
	regions = list(set(regions))

	# Logging frequency
	stepInterval = log['state'].iloc[1] - log['state'].iloc[0]

	# Dataframe with the root location
	logRootLocation = None
	if 'rootLocation' in log.columns.tolist():
		logRootLocation = log[['rootLocation']].reset_index(drop=True) 
		logRootLocation.rootLocation = logRootLocation.rootLocation.str.replace(' ', '') # Remove spaces if any
		log.drop(columns = 'rootLocation', inplace = True)

	# Add age root for basta
	if "basta" in runType.lower():
		log["ageRoot"] = getBastaAgeRoot(fileName, directory, runType, nSim, chainLength, stepInterval, burnIns)


	##############################################
	# For structured coalescent approximations
	# Transform backwards-in-time into forwards-in- 
	# time migration rates  
	##############################################
	if forwards.lower() in ["true", "all"] and "dta" not in runType.lower():
		log = forwardsRates(log, forwards, runType)


	##############################################
	# Summary statistics and ESS of all parameters
	##############################################
	# Summary statistics (min, mean, median, max, 95% CI)
	all_percentiles = [x for x in range(5, 100, 5)] + [1, 2.5, 99, 97.5]
	all_percentiles = sorted([x/100 for x in all_percentiles])
	results = log.drop(columns=['state']).describe(percentiles = all_percentiles).transpose(copy = True)
	results.drop(['count'], axis = 1, inplace = True)

	# ESS of all continuous parameters 
	# Migration rates are excluded when BSSVS is implemented
	discrete_or_rates = [x for x in log.columns if x in ["state", "sumNonZeroRates", "allTransitions"]] + \
						[x for x in log.columns if "indicators" in x or "nMigration" in x]

	if bssvs:
		discrete_or_rates += [x for x in log.columns if "backwards" in x or "forwards" in x]

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

	##############################################
	# Deal with spatial data when BSSVS is 
	# implemented and compute conditional ess  
	##############################################
	if bssvs:
		#bssvsResults, p_values = bssvsStatistics(
		bssvsResults = bssvsStatistics(
			log, 
			nRegions, 
			symmetric, 
			runType,
			stepInterval,  
			forwards
			)

		results = pd.concat([results, bssvsResults], sort = False, ignore_index = True) 

		# Remove lines corresponding to raw migration rates
		results = results[~results.parameter.str.match(r'^backwards_.*$|^forwards_.*$|^nMigration_.*$')]
		results.reset_index(inplace = True, drop = True)


	##############################################
	# Add root location if wanted
	##############################################
	if rootLocation:
		rootResults = rootLocationAnalysis(
			fileName, 
			directory, 
			runType, 
			regions, 
			logRootLocation, 
			regionDict,
			burnIn
			)

		results = pd.concat([results, rootResults], sort = False, ignore_index = True)

	##############################################
	# Add simulation number, sample size, protocol 
	# and matrix
	##############################################
	results.rename(columns={"2.5%":"2_5", "97.5%":"97_5"}, inplace = True)
	results.columns = results.columns.str.replace("%", "")
	results['nSim'] = nSim
	results['nSeq'] = nSeq
	results['protocol'] = protocol
	results['matrix'] = nMatrix
	results['model'] = runType

	return(results)

	










def getBastaRegionDict(fileName, directory):
	"""
	"""
	# Read comment lines of Basta Log file
	with open(directory + fileName, 'r') as commentLines:
		comments = takewhile(lambda s: s.startswith("#"), commentLines) 
		comments = list(comments)
	
	# Get regions in the correct order from the trait element 
	# It corresponds to the input data in the xml file
	for l in comments: 
		if '<trait id="typeTraitSet"' in l:
			regionList = re.findall(r"\.[0-9]*=(\w*),", l)

	regions = []
	for x in regionList:
		if x not in regions:
			regions.append(x)

	# Dictionary with keys as integers and values as region names
	out = dict(zip([str(x) for x in range(len(regions))], regions))

	return(out)










def renameColumns(log, runType, fileName, bastaRegionDict = None):
	"""
	Function that renames the columns and slice the beast log file according to
	the type of beast model.
	Arguments :
		- log (dataframe): log file corresponding to the beast .log.txt file without the burnin
		- runType (string): type of beast run
		For each runType, a specific dictionary is called from the beastLogDictionaries module
			"beast1" or "DTA" uses beast1Dict
			"mascot" used mascotDict
			"basta" uses bastaDict
		- fileName (str): name of the file

	It returns the same log file (values are not modified) with renamed columns.
	In the case of beast1 models, four columns are dropped. 
		- default.meanRate is the same as default.clock.rate for strict molecular clock models
		- regions.meanRate is the same as regions.clock.rate for strict molecular clock models
		- default.branchRates and regions.branchRates are equal to 0
	"""
	if 'dta' in runType.lower() :
		# Drop specific columns
		log.drop(columns = ['default.meanRate', 
			'regions.meanRate', 
			'default.branchRates',
			'regions.branchRates',
			'kappa.1', 
			'frequencies1.1', 
			'frequencies2.1', 
			'frequencies3.1', 
			'frequencies4.1'], 
			inplace = True,
			errors='ignore')

		# Rename columns
		log.rename(columns = beast1Dict, inplace = True)

	elif 'mascot' in runType.lower():
		# Add treeLikelihood key to mascotDict
		mascotDict['treeLikelihood.' + fileName.replace('.log.txt', '')] = "treeLikelihood"

		if any(["indicator" in x for x in log.columns]):
			# Add indicator keys to mascotDict
			ind = [x for x in log.columns if 'migration.indicator' in x]
			names = [x.replace('b_migration.', 'indicators_') for x in log.columns if 'b_migration' in x]
			names = [re.sub(r'_to_|_and_', '_', x)  for x in names]
			mascotDict.update(dict(zip(ind, names)))

		# Rename columns
		log.rename(columns = dict(mascotDict), inplace = True) 

	elif 'basta' in runType.lower():
		# Prepare column names for the structured coalescent variables
		columnNames = []
		newNames = []

		for col in log.columns:
			pattern = re.compile('|'.join(bastaDictElement.keys()))
			newName= pattern.sub(lambda x: bastaDictElement[x.group()], col)

			if 'nMigration' in newName:
				newName = newName.replace('to', '_')

			if 'likelihood' not in newName.lower() and 'hky' not in newName.lower():
				pattern = re.compile('|'.join(bastaRegionDict.keys()))
				newName = pattern.sub(lambda x: bastaRegionDict[x.group()], newName)

			if newName != col:
				newNames.append(newName)
				columnNames.append(col)

		# Append customed dictionary to static dictionary
		bastaDict.update( dict(zip(columnNames, newNames)) )

		# Rename columns
		log.rename(columns = bastaDict, inplace = True)

		# Replace integers by strings in rootLocation 
		log.rootLocation = log.rootLocation.astype('str').replace(bastaRegionDict)
		
		# Modify dtypes of nMigration columns
		col_str = [x for x in log.columns.tolist() if 'nMigration' in x]
		log[col_str] = log[col_str].replace('N', '-1')
		typeDict = dict(zip(col_str, ['int']*len(col_str)))
		log = log.astype(typeDict)
		log[col_str] = log[col_str].replace(-1, np.nan)

	else:
		raise ValueError('{0} is not a correct runType. \
			It can take the following values: beast1, DTA, mascot or basta'.format(runType))

	return(log)










def forwardsRates(log, forwards, beastModel):
	"""
	If a beast log file corresponds to a Mascot run, migration rates are backwards in time.
	To compare them with forwards in time migration rates, they can be converted into 
	forwards in time migration rates.

	forwards_{i,j} = backwards_{j,i} * Ne_{j} / Ne_{i}

	It returns the dataframe with supplementary migration rate columns with the suffix forwards_.

	For Markov jumps in Basta, the source and destination need to be inverted because they 
	correspond to backwards-in-time Markov jumps.

	"""
	
	out = log.copy()
	# allRates = [x for x in log.columns if "rates_" in x]

	for col in log.columns:
		if 'backwards_' in col:
			j = re.match(r'^backwards_(.*)_.*$', col).group(1)
			i = re.match(r'^.*_([^_]*)$', col).group(1)
			out['forwards_' + i + "_" + j] = log[col] * log['Ne_' + j] / log['Ne_' + i]
			
		if 'nMigration_' in col and "basta" in beastModel.lower():
			j = re.match(r'^nMigration_(.*)_.*$', col).group(1)
			i = re.match(r'^.*_([^_]*)$', col).group(1)
			out['nMigration_' + i + "_" + j] = log[col] 

	if forwards == "true":
		colsToDrop = [x for x in out.columns if 'backwards_' in x]
		out.drop(columns = colsToDrop, inplace=True)

	# if len(allRates) == 1:
	# 	out.drop(columns = allRates[0].replace("rates_", "nMigration"), inplace = True)

	return(out)










def bssvsStatistics(
	log, 
	nRegions, 
	symmetric, 
	runType, 
	stepInterval,
	forwards = "false", 
	essRates = pd.DataFrame()
	):
	"""
	If a beast log file corresponds to a BSSVS run, summary statistics (mean, max, min, quantiles, 
	95%-hpd interval, std) are computed on migration rates when the indicator variable is equal to 1.
	Bayes Factor (BF) supporting the migration rate is computed as well.

	Arguments:
		- log (pd.DataFrame): BEAST log file
		- nRegions (int): number of unique locations present in the log file
		- symmetric (bool): is the migration matrix symmetric or asymmetric 
		- runType (str): BEAST model 'dta', 'mascot_v*', 'basta'
		- forwards (str): argument for mascot beast files only
			if "true", then backwards-in-time migration rates are converted into forwards-in-time migration time 
			if "false", then migration rates are not converted 
			if "all", then both backwards-in-time and forwards-in-time migration rates are analyzed
		- essRates (pd.DataFrame): dataframe containing the ESS values of the raw migration rates
	
	Description of the output:
		- BF is computed using the bayesFactorRates function in the mathFunctions module. 
		It is directly based on Lemey et al., 2009.
	
		- Mean, min, max, quantiles and std are computed using the .describe() method from pandas.

		- 95%-HPD is computed using the hod function in the mathFunctions module.
		It is based on the hpd function from Beast math functions.
	"""

	##############################################
	# Determine parameters to compute
	# if forwards.lower() == "all":
	# 	prefix = ['rates']
	# 	if "dta" not in runType.lower():
	# 		prefix.append("forwards")
	
	# if forwards.lower() == "true":
	prefix = ['forwards']

	# elif forwards.lower() == "false":
	# 	prefix = ['rates']
	
	#if any('nMigration' == log.columns):
	if any(['nMigration_' in x for x in log.columns]):
		prefix.append('nMigration')

	##############################################
	# Indicator variables are coded separately or
	# in the migration parameter
	# It concerns only Mascot runs 
	separate = any(['indicators_' in x for x in log.columns]) 
		
	##############################################
	dataP = [] 
	kstest_dict = dict()

	for p in prefix: 
		cols = [x for x in log.columns if p + '_' in x]
		data = log[cols].copy()

		for colParameter in cols:
			
			if 'basta' in runType.lower() and p in ['forwards', 'nMigration']:
				source = re.match(r'^.*_([^_]*)$', colParameter).group(1)
				destination = re.match(r'^[a-zA-Z]*_(.*)_.*$', colParameter).group(1)
				colIndicator = "indicators_" + source + "_" + destination
			else:
				colIndicator = colParameter.replace(p, "indicators")

			# Summary statistics        
			if separate:
				posterior = log.loc[(log[colIndicator] == 1), colParameter]
			else :
				posterior = log.loc[(log[colParameter.replace("nMigration", "forwards")] > 0), colParameter]

			all_percentiles = [x for x in range(5, 100, 5)] + [1, 2.5, 99, 97.5]
			all_percentiles = sorted([x/100 for x in all_percentiles])
			posteriorDF = posterior.describe(percentiles = all_percentiles).to_frame().transpose().drop(columns=["count"])
			posteriorDF['2_5_hpd'] = hpd(list(posterior), lower_only = True) 
			posteriorDF['97_5_hpd'] = hpd(list(posterior), upper_only = True)

			# Convergence diagnostics
			if p == "nMigration" or len(posterior) == 0:
				posteriorDF['ESS'] = np.nan
			else :
				posteriorDF['ESS'] = ESS(list(posterior), stepInterval)

			# BF 
			if separate:
				posteriorDF['BF'] = bayesFactorRates(log[colIndicator], nRegions, sym = symmetric)
			else:
				x = mean(log[colParameter.replace("nMigration", "forwards")] > 0)
				if x == 1:
					x = x - (1 / log.shape[0])
				posteriorDF['BF'] = bayesFactorRates(x, nRegions, sym = symmetric)

			posteriorDF.reset_index(inplace = True, drop = True)

			# Change name of parameter 
			posteriorDF["parameter"] = 'bssvs_' + colParameter
			posteriorDF['p_value'] = np.nan
			if posteriorDF.BF.loc[0] >= 3 and posteriorDF.ESS.loc[0] < 200 :	
				# If ESS < 200 and BF >= 3
				# Perform a Kolmogorov Smirnov test between the start 
				# of the posterior and the end of the posterior to determine
				# if adding more samples from the posterior modifies the shape 
				# of the posterior
				# If not, we can consider that the parameter has mixed enough and 
				# it's OK if the ESS is low

				# KS test
				posterior.reset_index(inplace=True, drop=True)
				mid_index = math.ceil(posterior.shape[0] /2)
				x1 = posterior[range(mid_index)]
				posteriorDF["p_value"] = kstest(x1,posterior,alternative="two-sided").pvalue

			# Store subset dataframes
			dataP.append(posteriorDF)


	# Concatenate all dataframes
	out = pd.concat(dataP, axis = 0).reset_index(drop = True)

	return out










def posteriorDistnMig(log,runType):
	"""
	If a beast log file corresponds to a BSSVS run, summary statistics (mean, max, min, quantiles, 
	95%-hpd interval, std) are computed on migration rates when the indicator variable is equal to 1.
	Bayes Factor (BF) supporting the migration rate is computed as well.

	Arguments:
		- log (pd.DataFrame): BEAST log file
		- nRegions (int): number of unique locations present in the log file
		- symmetric (bool): is the migration matrix symmetric or asymmetric 
		- runType (str): BEAST model 'dta', 'mascot_v*', 'basta'
		- forwards (str): argument for mascot beast files only
			if "true", then backwards-in-time migration rates are converted into forwards-in-time migration time 
			if "false", then migration rates are not converted 
			if "all", then both backwards-in-time and forwards-in-time migration rates are analyzed
		- essRates (pd.DataFrame): dataframe containing the ESS values of the raw migration rates
	
	Description of the output:
		- BF is computed using the bayesFactorRates function in the mathFunctions module. 
		It is directly based on Lemey et al., 2009.
	
		- Mean, min, max, quantiles and std are computed using the .describe() method from pandas.

		- 95%-HPD is computed using the hod function in the mathFunctions module.
		It is based on the hpd function from Beast math functions.
	"""
	
	dataP = [] 

	p = 'nMigration' 
	separate = any(['indicators_' in x for x in log.columns]) 
	cols = [x for x in log.columns if p + '_' in x]
	data = log[cols].copy()

	for colParameter in cols:
		
		if 'basta' in runType.lower() :
			source = re.match(r'^.*_([^_]*)$', colParameter).group(1)
			destination = re.match(r'^[a-zA-Z]*_(.*)_.*$', colParameter).group(1)
			colIndicator = "indicators_" + source + "_" + destination
		else:
			colIndicator = colParameter.replace(p, "indicators")

		# Summary statistics        
		if separate:
			posterior = log.loc[(log[colIndicator] == 1), colParameter]
		else :
			posterior = log.loc[(log[colParameter.replace("nMigration", "forwards")] > 0), colParameter]

		posterior_df = posterior.to_frame(name = "value")
		posterior_df.reset_index(drop = True, inplace=True)
		posterior_df["parameter"] = colParameter

		# Store subset dataframes
		dataP.append(posterior_df)

	# Concatenate all dataframes
	out = pd.concat(dataP, axis = 0).reset_index(drop = True)
	out["model"] = runType

	return out











def rootLocationAnalysis(fileName, directory, model, regions, logFile, regionDict, burnIn):
	"""
	Function which computes the Kullback-Leibler divergence of the root location. 
	For DTA runs, the function uses the column 'regions' of the BEAST log file.
	For MASCOT runs, the function uses the corresponding tree posterior distribution.

	Arguments: 
		- fileName (str): name of the BEAST log file currently analyzed
		- directory (str): directory where to find the BEAST log file and the mcc tree
		- model (str): type of beast run (beast1, dta, mascot, basta)
		- regions (list): list of regions represented in the run under analysis 
		- logFile (pd.DataFrame): dataframe of the root location, it is passed only for 'dta' and 'beast1' runs
		- regionDict (dict): dictionary whose items are integers and values are the corresponding 
		location strings 

	It returns the KL divergence in a pd.DataFrame format. 
	"""
	# According to the model implemented, compute root location
	# probabilities and KL divergence
	if model.lower() in ["beast1"]:

		# Change fileName suffix
		fileName = re.sub(r"\..*", ".regions.states.log", fileName)

		# Read file
		file = pd.read_csv(directory + fileName, 
			comment = "#",
			header = 0, 
			sep = "\t").iloc[1000:]

		# Root location probababilites and KL
		out = file.groupby("regions").size().to_frame()
		out = out.reset_index().rename(columns={"regions":"parameter", 0:"value"})
		out.parameter = out.parameter.str.replace(' ', '')
		out.value /= file.shape[0]
		kl = KLRootPrediction(out.value, len(regions), True)


	elif "dta" in model.lower():

		# Root location probabilities
		out = logFile.groupby("rootLocation").size().to_frame()
		out = out.reset_index().rename(columns={"rootLocation":"parameter", 0:"value"})
		out.value = out.value / logFile.shape[0]

		# KL divergence 
		kl = KLRootPrediction(out.value, len(regions), True)


	elif "mascot" in model.lower():
		
		ann = annotationDict[model]
		rootLocation = []

		# Change filename suffix
		treeFiles = [directory + re.sub(r"\.log.*", ".trees.txt", fileName)]
		if os.path.exists(directory + re.sub(r"\.log.*", "_b.trees.txt", fileName)):
			treeFiles.append(directory + re.sub(r"\.log.*", "_b.trees.txt", fileName))

		# Get root location from tree posterior distribution
		treeCollectionYielder = dp.Tree.yield_from_files(
			files = treeFiles, 
			schema = 'nexus',  
			preserve_underscores = True,  
			extract_comment_metadata = True
		)

		for curr, tree in enumerate(treeCollectionYielder):
			if curr < burnIn:
				continue
			rootLocation.append(tree.seed_node.annotations.get_value(name=ann))

		out = pd.DataFrame(pd.Series(rootLocation).value_counts() / len(rootLocation) )
		out.reset_index(inplace=True)
		out.columns = ["parameter", "value"]

		# KL divergence
		kl = KLRootPrediction(out.value, len(regions), True)


	elif 'basta' in model.lower():

		# Root location probabilities
		out = logFile.groupby("rootLocation").size().to_frame()
		out = out.reset_index().rename(columns={"rootLocation":"parameter", 0:"value"})
		out.value = out.value / logFile.shape[0]

		# KL divergence 
		kl = KLRootPrediction(out.value, len(regions), True)

	

	# Add sampled locations that were not inferred as root location
	if out.shape[0] != len(regions):
		for x in regions:
			if x not in list(out.parameter):
				out = out.append(pd.DataFrame({'parameter':[x], 'value':[0.0]}), 
					ignore_index = True)

	# Add locations that were not sampled
	allRegions = list(regionDict.values())
	if out.shape[0] != len(allRegions):
		for x in allRegions:
			if x not in list(out.parameter):
				out = out.append(pd.DataFrame({'parameter':[x], 'value':[0.0]}), 
					ignore_index = True)

	# Add 'root_' prefix to parameter names 
	out.parameter = 'root_' + out.parameter.astype(str)

	# Add columns
	out = out.append(pd.DataFrame({"parameter": ["KL"], "value": [kl]}), ignore_index = True)

	return(out)










def ratesAndAdjustedBF(fileName, directory, beastModel = "dta", burnIn = 10):
	"""
	Function to compute the BF and adjusted BF of transition rates in DTA runs
		- fileName (str) : name of the log file with the .log.txt extension
		- directory (str): location of the original run
		- beastModel (str): name of the beast model. This name will be replaced by 
		beastModel_adjustedBF to load the associated run
		- burnIn (float): burnIn in percentage 

	The function returns a pandas dataframe with.

	"""
	###########################################
	# Get informations
	nSim = int(re.match(r'.*sim(\d*)_', fileName).group(1))
	nSeq = int(re.match(r'.*_(\d*)\.log.*', fileName).group(1))
	protocol = re.match(r'.*sim\d*_(.*)_\d*\.log.*', fileName).group(1)

	# Directory 
	directory = checkDirectory(directory)
	dir1 = directory
	dir2 = directory.replace(beastModel, beastModel + "_" + "adjustedBF")

	# Load log files without burnin (10% by default)
	log1 = pd.read_csv(dir1 + fileName, sep = "\t", comment = "#")
	burnIn1 = round(log1.shape[0] * burnIn / 100)
	log1 = log1.iloc[burnIn1:]
	log1.columns = log1.columns.str.replace(' ', '')

	log2 = pd.read_csv(dir2 + fileName, sep = "\t", comment = "#")
	burnIn2 = round(log2.shape[0] * burnIn / 100)
	log2 = log2.iloc[burnIn2:]
	log2.columns = log2.columns.str.replace(' ', '')
	
	###########################################
	# Output dataframe
	indicatorCols = [x for x in log1.columns if 'indicators' in x]
	outCols = ["parameter", "mean", "median", "97_5_hpd", "2_5_hpd", "97_5", "2_5", "BF", "adjustedBF"]
	out = pd.DataFrame(index = range(len(indicatorCols)), columns = outCols)
	out.parameter = [re.sub(r'.*indicators\.', '', x) for x in indicatorCols]

	# Number of regions
	regions = [re.sub('^.*\.', '', x) for x in log1.columns.tolist() if 'regions.rates.' in x]
	nRegions = len(set(regions))

	# Compute BF 
	for col in indicatorCols:

		if col not in log2.columns.tolist():
			raise ValueError("In condition {0}, column {1} is missing in the adjusted BF run".format(fileName, col))

		maskOut = out.parameter == re.sub(r'regions\.indicators\.', '', col)
		maskIn = log1[col] == 1
		migCol = col.replace('regions.indicators.', 'c_counts_')
		migCol = migCol.replace(".", "_") + "[1]"
		
		# Fill summary statistics of the transition rate
		out.loc[maskOut, "mean"] = log1[[migCol]][maskIn].mean()[0]
		out.loc[maskOut, "median"] = log1[[migCol]][maskIn].median()[0]
		out.loc[maskOut, "97_5"] = log1[[migCol]][maskIn].quantile(q = 0.975)[0]
		out.loc[maskOut, "97_5_hpd"] = log1[[migCol]][maskIn].apply(lambda x: hpd(x, upper_only = True))[0]
		out.loc[maskOut, "2_5"] = log1[[migCol]][maskIn].quantile(q = 0.025)[0]
		out.loc[maskOut, "2_5_hpd"] = log1[[migCol]][maskIn].apply(lambda x: hpd(x, lower_only = True))[0]

		# BF
		out.loc[maskOut, "BF"] = bayesFactorRates(log1[col].tolist(), nRegions)

		# Adjusted Bayes Factor
		p1 = log1[[col]].mean()[0]
		p2 = log2[[col]].mean()[0]
		out.loc[maskOut, "adjustedBF"] = (p1/(1-p1))/(p2/(1-p2))

	###########################################
	# Additionnal informations
	out.parameter = out.parameter.str.replace('.', '_')
	out['nSeq'] = nSeq
	out['nSim'] = nSim
	out['protocol'] = protocol

	return(out)









def getBastaAgeRoot(fileName, directory, runType, nSim, chainLength, stepInterval, burnIn):
	
	# Get mrst
	fasta_in = directory.replace(runType, "files") + fileName.replace("sim" + str(nSim) + "_", "sim" + str(nSim) + "_sequences_").replace(".log.txt", ".fasta")
	record = SeqIO.index(fasta_in, 'fasta')
	taxa = list()
	with open(fasta_in) as in_handle:
		for title, seq in SimpleFastaParser(in_handle):
			taxa.append(title)

	dates = [float(re.search('(?<=_)[0-9]+[.]*[0-9]+$', t).group(0)) for t in taxa]
	mrst = max(dates)

	# Get age root
	allAgeRoot = []
	treeFiles = [directory+f for f in os.listdir(directory) if re.match(fileName.replace(".log.txt", "")+".*trees.txt", f)]

	for i, f in enumerate(treeFiles):
		treeLabels = []
		ageRoot = []
		post_burn_in_labels = range(0, chainLength[i] * stepInterval, stepInterval)[burnIn[i]: chainLength[i]]

		treeCollectionYielder = dp.Tree.yield_from_files(
			files = [f], 
			schema = 'nexus',  
			preserve_underscores = True,  
			extract_comment_metadata = True
			)

		for tree in treeCollectionYielder:
			label = float(tree.label.replace("STATE_", ""))
			if label in post_burn_in_labels:
				treeLabels.append(label)
				ageRoot.append(mrst - max(tree.calc_node_ages(is_force_max_age =True)))

		# Add nan when some trees are missing
		for ind in [j for j, e in enumerate(post_burn_in_labels) if e not in treeLabels]:
			ageRoot.insert(ind, np.nan)

		allAgeRoot += ageRoot

	return(allAgeRoot)

