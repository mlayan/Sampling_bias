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
import pandas as pd
import numpy as np
import dendropy

## Import custom functions
from mathFunctions import *
from simulatedTrees import getAncestors
from handleDirectories import checkDirectory
from itertools import takewhile

## Import regions and log files columns dictionaries
from beastLogDictionaries import * 







def logFileWrangler(fileName, directory, nMatrix, runType, regionDict,
	rootLocation = False, forwards = "false", ess = False):
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
	log = pd.read_csv(directory + fileName, sep = "\t", comment = "#").iloc[1000:]

	# Remove spaces and rename columns
	log.columns = log.columns.str.replace(' ', '')

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
	nSim = int(re.match(r'sim(\d*)_', fileName).group(1))
	nSeq = int(re.match(r'.*_(\d*)\.log.*', fileName).group(1))
	protocol = re.match(r'sim\d*_(.*)_\d*\.log.*', fileName).group(1)
	
	# Number of regions 
	regions = [re.sub('^.*_', '', x) for x in log.columns.tolist() \
			if 'rates_' in x]
	nRegions = len(set(regions))

	# Symmetric or asymmetric matrix 
	if len(regions) == (nRegions * (nRegions-1)) / 2:
		symmetric = True
	elif len(regions) == nRegions * (nRegions - 1):
		symmetric = False
	elif len(regions) == 1 and nRegions == 1:
		symmetric = True
		regions = regions = [x for x in log.columns.tolist() if 'rates_' in x][0]
		regions = regions.split("_")
		regions.remove("rates")
		nRegions = len(regions)
	else:
		raise ValueError("The number of migration rates ({0}) doesn't correspond to \
							the number of regions ({1})".format(len(regions), nRegions))

	# List of regions
	regions = list(set(regions))

	# Dataframe with the root location
	logRootLocation = None
	if 'rootLocation' in log.columns.tolist():
		logRootLocation = log[['rootLocation']].reset_index(drop=True) 
		logRootLocation.rootLocation = logRootLocation.rootLocation.str.replace(' ', '') # Remove spaces if any
		log.drop(columns = 'rootLocation', inplace = True)

	##############################################
	# Mascot runs only
	# Transform backwards in time in forwards in 
	# time migration rates  
	##############################################
	if forwards.lower() in ["true", "all"]:
		if "mascot" not in runType.lower() and 'basta' not in runType.lower():
			raise ValueError('Forwards in time migration rates can be computed only for Mascot runs. \
				The "forwards" argument has not a correct value for a {0} model'.format(runType))
		else:
			log = forwardsRates(log)

	if forwards.lower() == "true":
		cols = [c for c in log.columns if 'rates_' in c]
		log.drop(columns=cols, inplace=True)
		

	##############################################
	# Get describe(), the 95% HPD intervals and 
	# ESS from all parameters
	##############################################
	# Describe() on log
	results = log.drop(columns=['state']).describe(percentiles = 
		[0.025, 0.5, 0.975]).transpose(copy = True)
	results.drop(['count'], axis = 1, inplace = True)

	# ESS() on log
	if ess:
		stepInterval = log['state'].iloc[1] - log['state'].iloc[0]
		results['ESS'] = log.drop(columns=['state']).apply(
			lambda x: ESS(x, stepInterval)
			)

	# Compute 95% HPD intervals 
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
	# implemented
	##############################################
	# Is BSSVS implemented?
	for x in log.columns.tolist():

		if 'indicators' in x:
			if ess:
				# essRates = results[['ESS', 'parameter']].copy()
				# essRates.parameter.replace('rates_', 'backwards_', regex=True, inplace=True)

				essRates = results.loc[results.parameter.str.match('^rates_|^forwards_'), ['ESS', 'parameter']]            
				if 'dta' in runType.lower():
					essRates.parameter.replace('rates_', 'forwards_', regex=True, inplace=True)
				elif 'mascot' in runType.lower(): 
					essRates.parameter.replace('rates_', 'backwards_', regex=True, inplace=True)

				essRates['parameter'] = 'bssvs_' + essRates['parameter'].astype(str)
				bssvsResults = bssvsStatistics(log, nRegions, symmetric, runType, forwards, essRates)

			else:
				bssvsResults = bssvsStatistics(log, nRegions, symmetric, runType, forwards)

			
			# Add BF to Markov Jumps 
			results['param'] = results.parameter.str.replace(r'^[a-zA-Z]*_', '')
			results.loc[~results.parameter.str.contains('nMigration'), 'param'] = np.nan
			results = pd.merge(results, bssvsResults[['BF', 'param']], on = ['param'], how = "left")

			# Concatenate results 
			results = pd.concat([results, bssvsResults], sort = False, ignore_index = True) 
			results.drop(columns = "param", inplace = True)

			break

	##############################################
	# Add root location if wanted
	##############################################
	if rootLocation:
		rootResults = rootLocationAnalysis(fileName, directory, runType, regions, 
			logRootLocation, regionDict)

		results = pd.concat([results, rootResults], sort = False, ignore_index = True)

	##############################################
	# Add simulation number, sample size, protocol 
	# and matrix
	##############################################
	results.rename(columns={"2.5%":"2_5", "97.5%":"97_5", "50%":"50"}, inplace = True)
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
			if 'likelihood' not in newName:
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
		log[col_str] = log[col_str].replace('N', '0')
		typeDict = dict(zip(col_str, ['int']*len(col_str)))
		log = log.astype(typeDict)

	else:
		raise ValueError('{0} is not a correct runType. \
			It can take the following values: beast1, DTA, mascot or basta'.format(runType))

	return(log)










def forwardsRates(log):
	"""
	If a beast log file corresponds to a Mascot run, migration rates are backwards in time.
	To compare them with forwards in time migration rates, they can be converted into 
	forwards in time migration rates.

	forwards_{i,j} = backwards_{j,i} * Ne_{j} / Ne_{i}

	It returns the dataframe with supplementary columns with the suffix forwards_.
	"""
	
	for col in log.columns:
		if 'rates_' in col:
			j = re.match(r'^rates_(.*)_.*$', col).group(1)
			i = re.match(r'^.*_([^_]*)$', col).group(1)
			log['forwards_' + i + "_" + j] = log[col] * log['Ne_' + j] / log['Ne_' + i]

	return(log)










def bssvsStatistics(log, nRegions, symmetric, runType, forwards = "false", 
	essRates = pd.DataFrame()):
	"""
	If a beast log file corresponds to a BSSVS run, statistics (mean, max, min, quantiles, 
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
	# Compute BF based on indicator variable 
	# of migration rates
	colsI = [x for x in log.columns if 'indicators_' in x]
	bfD = log[colsI].apply(
		lambda x: bayesFactorRates(x, nRegions, sym = symmetric), 
		axis = 0
		)
	bfD = pd.DataFrame(bfD, columns = ['BF']).rename_axis('param').reset_index()
	bfD.param = bfD.param.str.replace('indicators_', '')
	
	##############################################
	# rates*indicators dataframe
	# Convert the dataframe in a long format with 
	# one columns with indicators and the other with rates
	if forwards.lower() == "all":
		merged = longFormat(log, ['forwards', 'rates'])
	
	elif forwards.lower() == "true":
		merged = longFormat(log, 'forwards')

	elif forwards.lower() == "false" and "dta" in runType.lower():
		merged = longFormat(log, 'rates')

	# Filter the dataframe on the indicator to keep only 
	# row with indicators equal to 1
	merged = merged[merged['indicators'] == 1].\
			drop(columns = ['id', 'indicators']).reset_index(drop = True).copy()

	# Dataframe of the summary statistics for spatial rates
	bssvsStats = merged.groupby(['param', 'parameter']).describe(percentiles = [0.025, 0.5, 0.975])
	bssvsStats.columns = bssvsStats.columns.get_level_values(1)

	# Compute the 95% HPD 
	bssvsHpd = merged.groupby(['param', 'parameter'], as_index = False).\
	agg([lambda x : hpd(x, upper_only = True), lambda x: hpd(x, lower_only = True)])

	bssvsHpd.columns = bssvsHpd.columns.get_level_values(1)
	bssvsHpd.columns = ['97_5_hpd', '2_5_hpd']

	# Concatenate hpd and stats dataframes
	bssvs = pd.concat([bssvsHpd, bssvsStats], axis = 1, sort = False)
	bssvs.drop(columns = 'count', inplace = True) # Remove count column
	bssvs.reset_index(inplace = True) # convert row indexes into columns

	# Concatenate dataframes
	out = pd.merge(bfD, bssvs, on = ['param'])
	#out.drop(columns = 'parameter', inplace = True)
	#out.rename(columns = {'param':'parameter'}, inplace = True)

	# Add the ESS if they were calculated 
	if not essRates.empty:
		out = pd.merge(out, essRates, on = ['parameter'])

	del merged, bssvsHpd, bssvsStats, bfD, bssvs

	return(out)










def longFormat(log, prefix):
	"""
	Part of the bssvsStatistics function which converts the log dataframe (wide format) into a long format.
		- log (pandas dataframe) : beast log file with standardized column names
		- prefix (list of strings) : prefix of the columns to be multiplied with its corresponding 
		indicator variable

	It returns the dataframe in long format.
	"""
	#prefix = ['forwards', 'rates']
	if type(prefix) == str:
		prefix = [prefix]

	dataP = [] # List of dataframes to concatenate

	for p in prefix:
		# Subset the dataframe to columns containing the prefix
		cols = [x for x in log.columns if p in x]
		data = log[cols].copy()
		data['id'] = data.index

		# Rearrange the dataframe in long format
		data = data.melt(id_vars = ['id'], var_name = 'parameter', value_name = p)

		if p == "rates":
			if len(prefix) > 1:
				data.parameter.replace(to_replace = p, 
					value = 'bssvs_backwards', 
					inplace = True, 
					regex = True)
			else:
				data.parameter.replace(to_replace = p, 
					value = 'bssvs_forwards', 
					inplace = True, 
					regex = True)

			data['param'] = data.parameter.replace(to_replace='^[a-z\_]*', 
				value = '', regex = True) 
			dataP.append(data)

		else:
			data.parameter.replace(to_replace = p, 
				value = 'bssvs_forwards', 
				inplace = True, 
				regex = True)
			data.rename(columns = {'forwards':'rates'}, inplace = True)
			data['destination'] = data.parameter.replace(to_replace='^\w*_', 
														value = '', regex = True)
			data['source'] = data.parameter.replace(to_replace='^[a-z\_]*|_\w*$', 
														value = '', regex = True)
			data['param'] = data.destination + "_" + data.source
			data.drop(columns=['destination', 'source'], inplace=True)
			dataP.append(data)


	dataP = pd.concat(dataP, sort = False, ignore_index = True)

	# Do the same with the indicator dataframe
	colsI = [x for x in log.columns if 'indicators_' in x]
	dataI = log[colsI].copy()
	dataI['id'] = dataI.index
	dataI = dataI.melt(id_vars = ['id'], var_name = 'param', value_name = 'indicators')
	dataI.param = dataI.param.replace(to_replace = '^[a-z]*_', value = '', regex = True)

	# Merge dataframes
	merged = pd.merge(dataI, dataP, on =  ['id', 'param'])

	return(merged)










def rootLocationAnalysis(fileName, directory, model, regions, 
	logFile, regionDict):
	"""
	Function which computes the Kullback-Leibler divergence of the root location. 
	For DTA runs, the function uses the column 'regions' of the BEAST log file.
	For MASCOT runs, the function uses the corresponding mcc tree file.

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
		# Change fileName suffix
		fileName = re.sub(r"\.log.*", ".mcc.tree", fileName)

		# Read nexus tree
		tree = dendropy.Tree.get(path = directory + fileName, 
			schema = 'nexus', 
			preserve_underscores = True)

		# Root annotations
		rootAnnotations = tree.seed_node.annotations
		# nRegions = len(rootAnnotations.get_value(name="max.set"))
		# if nRegions != len(regions):
		# 	nRegions = len(regions)

		# KL divergence
		p = [rootAnnotations.get_value(name = x) for x in regions]
		p = [float(x) if x else 0.0 for x in p]
		# p = [float(x) for x in rootAnnotations.get_value(name = "max.set.prob")]
		kl = KLRootPrediction(p, len(regions), True)

		# Probability per region 
		# stored in a dataframe
		# r = [re.sub('"', '', x) for x in rootAnnotations.get_value(name = "max.set")]
		# r = [re.sub(' ', '', x) for x in r]

		# # Add regions in the sample but not inferred as root location
		# for x in regions:
		# 	if x not in r:
		# 		r.append(x)
		# 		p.append(0.0)

		# # Add regions not in the sample with NA as value
		# for x in allRegions:
		# 	if x not in r:
		# 		r.append(x)
		# 		p.append(float('nan'))

		# Create the dataframe
		out = pd.DataFrame({"parameter": regions, "value": p})


	elif 'basta' in model.lower():
		# if not bastaRegionDict:
		# 	raise ValueError("A dictionary for regions needs to be passed for basta outputs")

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










def ratesAndAdjustedBF(fileName, directory, beastModel):
	"""
	Function to compute the BF and adjusted BF of transition rates in DTA runs
		- fileName (str) : name of the log file with the .log.txt extension
		- directory (str): location of the original run
		- beastModel (str): name of the beast model. This name will be replaced by 
		beastModel_adjustedBF to load the associated run

	The function returns a pandas dataframe with.

	"""
	###########################################
	# Get informations
	nSim = re.sub(r'^sim|_.*$', '', fileName)
	nSeq = re.sub(r'^.*_|\.log\.txt', '', fileName)
	protocol = re.sub('^sim[0-9]*_|_[0-9]*\.log\.txt$', '', fileName)

	# Directory 
	directory = checkDirectory(directory)
	dir1 = directory
	dir2 = directory.replace(beastModel, beastModel + "_" + "adjustedBF")

	# Load log files without burnin (considered 10%)
	log1 = pd.read_csv(dir1 + fileName, sep = "\t", comment = "#").iloc[1000:]
	log1.columns = log1.columns.str.replace(' ', '')
	log2 = pd.read_csv(dir2 + fileName, sep = "\t", comment = "#").iloc[1000:]
	log2.columns = log2.columns.str.replace(' ', '')
	
	###########################################
	# Output dataframe
	indicatorCols = [x for x in log1.columns if 'indicators' in x]
	outCols = ["parameter", "mean", "median", "97_5_hpd", "2_5_hpd", "97_5", "2_5", "BF", "adjustedBF"]
	out = pd.DataFrame(index = range(len(indicatorCols)), columns = outCols)
	out.parameter = [re.sub(r'.*indicators\.', '', x) for x in indicatorCols]

	# Number of regions
	regions = [re.sub('^.*\.', '', x) for x in log1.columns.tolist() if 'rates' in x]
	nRegions = len(set(regions))

	# Compute BF 
	for col in indicatorCols:

		maskOut = out.parameter == re.sub(r'.*indicators\.', '', col)
		maskIn = log1[col] == 1
		rateCol = col.replace('indicators', 'rates')
		
		# Fill summary statistics of the transition rate
		out.loc[maskOut, "mean"] = log1[[rateCol]][maskIn].mean()[0]
		out.loc[maskOut, "median"] = log1[[rateCol]][maskIn].median()[0]
		out.loc[maskOut, "97_5"] = log1[[rateCol]][maskIn].quantile(q = 0.975)[0]
		out.loc[maskOut, "97_5_hpd"] = log1[[rateCol]][maskIn].apply(lambda x: hpd(x, upper_only = True))[0]
		out.loc[maskOut, "2_5"] = log1[[rateCol]][maskIn].quantile(q = 0.025)[0]
		out.loc[maskOut, "2_5_hpd"] = log1[[rateCol]][maskIn].apply(lambda x: hpd(x, lower_only = True))[0]

		# BF
		out.loc[maskOut, "BF"] = log1[[col]].apply(lambda x: bayesFactor(x, nRegions))[0]

		# Adjusted Bayes Factor
		p1 = log1[[col]].sum()[0] / log1.shape[0]
		p2 = log2[[col]].sum()[0] / log2.shape[0]
		out.loc[maskOut, "adjustedBF"] = (p1/(1-p1))/(p2/(1-p2))


	###########################################
	# Additionnal informations
	out.parameter = out.parameter.str.replace('.', '_')
	out['nSeq'] = nSeq
	out['nSim'] = nSim
	out['protocol'] = protocol

	return(out)