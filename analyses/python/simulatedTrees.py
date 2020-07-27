#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MODULE TO ANALYZE SIMULATED TREES
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-02-27' 
__last_update__ = '2020-04-21'

# Modules
import os
import sys
import re
import random
import dendropy
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from Bio import Phylo
from io import StringIO
from statistics import mean

# Import fixed data from trees
from simulatedTreesData import *
from handleDirectories import checkDirectory








 
def getAncestors(individual, transmissionChain):
	"""Function to get the ancestral lineage of individual from a simulated epidemic 
		- individual: name of the individual for which we want to get the ancetsral lineage
		- transmissionChain: dataframe corresponding to the transmission chain. 
		It should contain two columns:
			- case: name of the cases in the following format "patch_uniqueid"
			- source: name of the infector of case in the same format "patch_uniqueid"

	This function returns the ancestral lineage in a list
	"""

	##################################################
	## Get all the ancestors of the individual
	ancestorsList = []
	ancestor = None

	while not ancestor in root:
		ancestor = transmissionChain.loc[transmissionChain['case'] == individual, 'source'].values[0]
		ancestorsList.append(ancestor)
		individual = ancestor

	return(ancestorsList)










def countsAndRates(df, t):
	""" Function to compute the number of migration events from a transmission chain dataframe
	Parameters:
		- df (pandas dataframe): transmission chain with columns
			regions : names of cases' region
			regions.source : names of sources' region
		- time (int) : time span of the transmission chain in days

	It returns a simplified dataframe with the following columns:
		- value (float) : number or rate of migration event  
		- parameterChar (str) :  
		- parameter (str) : name of the parameter, either migrationCount or migrationRate
	"""
	# Slice the dataframe
	df = df[['regions.source', 'regions']]

	# Group columns to get counts
	counts = df.groupby(['regions.source', 'regions']).size().to_frame(name = 'value').reset_index()

	# Transform source and region columns into parameterChar column
	counts['parameterChar'] = counts['regions.source'] + "_" + counts['regions']
	counts['parameterChar'] = counts['parameterChar'].str.replace(' ', '')
	counts.drop(columns = ['regions.source', 'regions'], inplace=True)

	# Add Na for migrations which do not occur in the simulation
	absent = []
	absentValue = []
	for x in migrationList:
		if x not in list(counts.parameterChar):
			absent.append(x)
			absentValue.append(float("nan"))

	counts = counts.append(pd.DataFrame({'parameterChar': absent,
										'value': absentValue}), 
							ignore_index = True)

	# Create a copy to compute migration rates 
	rates = counts.copy()
	rates['value'] = rates['value'] / t

	# Add parameter column
	counts['parameter'] = "migrationCount"
	rates['parameter'] = "migrationRate"

	# Regions frequencies
	frequencies = df.groupby(['regions']).size().to_frame(name = 'value') / df.shape[0]
	frequencies = frequencies.reset_index().rename(columns = {'regions':'parameterChar'})
	frequencies.parameterChar = frequencies.parameterChar.str.replace(' ', '')
	frequencies['parameter'] = 'locationFrequency'

	# Concatenate dataframes
	out = pd.concat([counts, rates, frequencies], ignore_index=True, sort = False)

	return(out)

	"""
	# Transform source and region columns into parameterChar column
	counts['regions.source']= counts['regions.source'].str.replace(' ', '')
	counts['regions']= counts['regions'].str.replace(' ', '')
	counts['parameterChar'] = counts['regions.source'] + "_" + counts['regions']
	counts.drop(columns = ['regions.source', 'regions'], inplace=True)

	# Create two copies for the rates 
	rates = counts.copy()
	rates['value'] = rates['value'] / t

	# Add parameter column
	counts['parameter'] = "migrationCount"
	rates['parameter'] = "migrationRate"

	# Concatenate dataframes
	out = pd.concat([counts, rates], ignore_index=True)

	return(out)
	"""









"""
def trimPhylogeny(subset, phylogeny, type):
	Function to trim a phylogeny based on a set of samples and the complete phylogeny.
		- subset: subset of samples to which the phylogeny should be trimmed
		- phylogeny: dataframe corresponding to the transmission chain.
		It should contain two columns:
			- case: name of the cases in the following format "patch_uniqueid"
			- source: name of the infector of case in the same format "patch_uniqueid"
			- ??
		- type: string specifying how the phylogeny should be trimmed
		It can take the following arguments:
			- date: returns the phylogeny with no sequences older than the oldest sequence in from subset
			- mrca: returns the phylogeny with the complete offspring of the MRCA of subset until the latest sequence in subset 
			- stringent: returns the phylogeny with the offspring of the MRCA restricted to branches ending by one of the sequence from subset

	It returns the phylogeny dataframe with rows corresponding to the desired ancestors
	"""










def treeScanNoLabel(ind, df):
	"""
	Recursive function that creates from a dataframe and a tip/node of a phylogenetic tree
	the subtree in newick format from the tip/node as the root.

	Arguments:
		- ind (string): name of the node or individual 
		- df (dataframe): phylogenetic tree in a dataframe format
			- column 1: name of nodes
			- column 2: name of nodes' children

	The function returns:
		- out (string): the newick tree 
	"""
	if len(df.columns) != 2:
		raise ValueError("Dataframe has not the correct dimension : {} columns instead of 2 or 3".format(len(df.columns)))

	# Rename the first two columns
	df.rename(columns = {list(df)[0]: 'source', list(df[1]):'child'}, inplace = True)           

	descCond = (df['source'] == ind)    # Condition to get the offspring of ind
	desc = list(df[descCond]['child']) # List of ind's offspring
	descN = sum(descCond)        # Size of ind's offspring

	out = ""

	# ind has at least one descendant
	if descN != 0:
		if descN == 1: # One descendant
			out += treeScan(desc[0], df, innerlabel)
			return "(" + out + "),"

		else: # Several descendants
			out2 = []
			for d in desc:
				innerDesc = treeScan(d, df, innerlabel)
				out2.append(innerDesc)
			return "(" + ",".join(out2) + ")"

	# ind is a tip of the tree
	else:
		out += ind
		return(out)










def treeScan(ind, df, nodeDate = None):
	"""
	Recursive function that creates from a dataframe and a tip/node of a phylogenetic tree
	the subtree in newick format from the tip/node as the root. Branch lengths are informed
	in the newick tree.

	Arguments:
		- ind (string): name of the node or individual 
		- df (dataframe): phylogenetic tree in a dataframe format
			- column 1: name of nodes
			- column 2: name of nodes' children
		- nodeDate (int) : date of the node if it is not present as 'child' in df

	The function returns:
		- out (string): the newick tree 
		- nodeDate (int): date of ind 
	"""

	if len(df.columns) != 3:
		raise ValueError("Dataframe has not the correct dimension : {} columns instead of 3".format(len(df.columns)))
	df.columns = ['source', 'child', 'label']

	if not nodeDate : 
		nodeDate = df[df['child'] == ind]['label'].values
		if len(nodeDate) > 1:
			raise ValueError("{0} appears multiple times in the dataframe as a tip".format(ind))
		nodeDate = nodeDate[0]


	descCond = (df['source'] == ind)    # Condition to get the offspring of ind
	desc = list(df[descCond]['child'])  # List of ind's offspring
	descN = sum(descCond)               # Size of ind's offspring

	# ind has at least one descendant
	if descN != 0:

		if descN == 1: # One descendant
			return treeScan(desc[0], df)

		else: # Several descendants
			out = []
			for d in desc:
				inner, innerLabel = treeScan(d, df)
				inner += ":" + str(innerLabel - nodeDate)
				out.append(inner)

			return "(" + ",".join(out) + ")", nodeDate

	# ind is a tip of the tree
	else:
		label = df[df['child'] == ind]['label'].values
		if len(label) > 1:
			raise ValueError("{0} appears multiple times in the dataframe as a tip".format(ind))

		if label[0].is_integer():
			labelString = round(label[0])
		else:
			labelString = round(label[0], 6)

		out = ind + "_" + str(labelString)
		return(out, label[0])










def newickTree(rootTree, df, rootSDate, innerlabel = True):
	"""
	Function that returns the final Newick tree from a phylogenetic tree in a dataframe format.
	It is based on the treeScan(NoLabel) functions.
	Arguments:
		- rootTree (string): name of the tree root
		- df (dataframe) : phylogenetic tree
			- column 1: nodes name
			- column 2: children name
			- column 3: sampling date of children
		- rootSDate (int): sampling date of rootTree
		- innerlabel (bool): add branch length in the Newick tree

	Output:
		- out (string): Final Newick tree
	"""
	if innerlabel:
		(out, nodeDate) = treeScan(ind = rootTree, df = df, nodeDate = rootSDate) 

		if nodeDate != rootSDate:
			out = "(" + out 
			out += ":" + str(nodeDate - rootSDate) + ");"
		else:
			out += ";"

		return(out)

	else:
		out = treeScanNoLabel(rootTree, df)
		out += ";"
		return(out)










def trimPhylogeny(transChain, tips, ancestorsList, ancestorsCommon, rootTree):
	"""
	Function that returns the subtree corresponding to the tree with tips as tips from a 
	larger phylogenetic tree in transChain.
	Arguments:
		- transChain (dataframe): phylogenetic tree
		- tips (list of strings): list of tip names
		- ancestorsList (list of lists of strings): list of ancestors names for each tip
		- ancestorsCommon (set of strings): set of ancestors common to all tips
		- rootTree (string): MCRA of all tips

	- Output:
		- out (dataframe): trimmed phylogeny

	"""
	## Get ancestors of the subtree 
	ancestorsAll = set().union(*ancestorsList)
	ancestorsSubtree = list(ancestorsAll - ancestorsCommon)
	ancestorsSubtree.append(rootTree)

	## Subset the transmission chain to ancestorsSubtree et subtree tips
	out = transChain[['source', 'case', 'samp.dates']]
	out = out[out.source.isin(ancestorsSubtree)]
	out = out[out.case.isin(ancestorsSubtree + tips)]

	# Additionnal nodes for sample individuals which are also part 
	## of the inner phylogeny
	for x in tips:
		if x in ancestorsSubtree:
			addNode = "a_" + x
			nodeDate = out[out['case'] == x ]['samp.dates'].iloc[0]        
			out.replace({x:addNode}, inplace = True)
			out = out.append({'case': x, 'source': addNode, 'samp.dates': nodeDate + 0.001}, ignore_index = True)

	return(out)










def trimPhylogeny2(transChain, tips, ancestorsList, ancestorsCommon, rootTree):
	"""
	Function that returns the subtree corresponding to the tree with tips as tips from a 
	larger phylogenetic tree in transChain.
	Arguments:
		- transChain (dataframe): phylogenetic tree
		- tips (list of strings): list of tip names
		- ancestorsList (list of lists of strings): list of ancestors names for each tip
		- ancestorsCommon (set of strings): set of ancestors common to all tips
		- rootTree (string): MCRA of all tips

	- Output:
		- out (dataframe): trimmed phylogeny

	"""
	## Get ancestors of the subtree 
	ancestorsAll = set().union(*ancestorsList)
	ancestorsSubtree = list(ancestorsAll - ancestorsCommon)
	ancestorsSubtree.append(rootTree)

	## Subset the transmission chain to ancestorsSubtree et subtree tips
	out = transChain[['source', 'case', 'samp.dates']]
	out = out[out.source.isin(ancestorsSubtree)]
	out = out[out.case.isin(ancestorsSubtree + tips)]

	# Additionnal nodes for sample individuals which are also part 
	## of the inner phylogeny
	for x in tips:
		if x in ancestorsSubtree:
			addNode = "a_" + x    # Additionnal node which will be considered as the ancestor of x
			nodeDate = out[out['case'] == x ]['samp.dates'].iloc[0]   # Sampling date of x     
			out.at[out['case'] == x, 'samp.dates'] = nodeDate - 0.000001
			out.replace({x:addNode}, inplace = True)   # All occurences of x are replaced by the additionnal node
			out = out.append({'case': x, 'source': addNode, 'samp.dates': nodeDate}, 
			ignore_index = True) # new row corresponding to the branch (additionnalnode -> x)

	return(out)










def migrationEvents(filename, matrix, protocol, 
	extractNewickTree = False, directory = None):
	"""
	Function that computes the number of migration events between locations on the complete 
	phylogeny or on the trimmed phylogeny from simulated transmission chains.

	Arguments:
		- filename (str): name of the simulated transmission chain
		- matrix (int): number of the spatial coupling matrix used to simulate the transmission chain
		- protocol (str): name of the sampling protocol 
			- "all": considers the complete phylogeny and doesn't write a Newick tree
			- see all protocols in the "protocols" variable of the simulatedTreesData module
		- extractNewickTree (bool): extract the newick phylogenetic tree from the dataframe and 
		write it to a .nwk file
		- directory (str): name of the directory where to write the Newick tree file  


	Outputs: 
		- It writes the Newick tree when a specific sampling protocol is passed.
		The file name is filename.nwk 

		- The dataframe of the number of migration events between two locations. 
	"""
	print(filename, protocol)

	##################################################  
	if extractNewickTree and protocol == "all":
		raise ValueError('No tree is extracted when migration events \
			are evaluated over the entire phylogeny.')

	##################################################
	## Directory name 
	if directory:
		if not re.search(r'/$', directory):
			directory += "/"
	else:
		directory = ""

	##################################################
	## Load and rearrange the transmission chain 
	# Load
	transChain = pd.read_csv(directory + filename, sep="\t", header = 0)

	# Drop the row corresponding to the root
	transChain = transChain[transChain['case'] != "0_0"]

	# Add regions.source column to transChain
	transChain['regions.source'] = transChain['patch.source'].replace(regionDic)

	# nSim 
	nSim = re.findall(r'sim([0-9]+)_', filename)[0]

	##################################################
	## Extract informations 
	if protocol == "all":
		# Number of days corresponding to the entire epidemic
		nDays = round(365.25*nYears)

		# nSeq doesn't take the size of the sample but na
		nSeq = float("nan")

		## Slice the transChain dataframe to get 
		##    regions and regions.source only
		t = transChain

	else:
		if protocol not in protocols:
			raise ValueError("The value of the protocol argument is not referenced. \
				\nReferenced values are 'all' and values in the simulatedTreesData module.")

		##################################################
		## Get the list of all the common ancestors of 
		## individuals that were sampled
		# List of individuals in the sample
		sample = transChain['case'][transChain[protocol] == 1]
		sample = list(sample)

		## Sampling date of the last sample
		lastDate = transChain['samp.dates'][transChain[protocol] == 1].sort_values(ascending = False).iloc[0]
		
		## Get all the ancestors 
		ancestors = [getAncestors(x, transChain[['case', 'source']]) for x in sample]

		## Intersection of ancestors 
		inter = set.intersection(*map(set, ancestors))
		interL = list(inter)
		if len(interL) == 1:
			if root == interL:
				rootS = root[0]
				rootSLocation = rootLocation
				rootSDate = rootDate
			else:
				raise ValueError("A single MCRA has been found : {0}. \
					It does not correspond to {1}.".format(interL, root))
		else:
			## Identification of sampling date and location of MRCA
			sub = transChain[transChain['case'].isin(interL)].sort_values(by='samp.dates', ascending = False)
			rootS = sub['case'].iloc[0]
			rootSLocation = sub['regions'].iloc[0]
			rootSDate = sub['samp.dates'].iloc[0]

		# nSeq
		nSeq = re.match(r'.*_(\d+)$', protocol).group(1)
		if int(nSeq) not in sampleSize:
			raise ValueError("{0} is not correctly formated.".format(protocol))

		# EpidemicDuration in days  
		treeHeight = lastDate - rootSDate
		nDays = round(treeHeight*365.35)

		## Trim phylogeny
		t = transChain[(transChain['samp.dates'] >= rootSDate) & (transChain['samp.dates'] <= lastDate)].copy().reset_index()

		##################################################
		if extractNewickTree:
			## Create the Newick tree
			df = trimPhylogeny2(t, sample, ancestors, inter, rootS)
			tree = newickTree(rootS, df, rootSDate)

			# Newick tree file name
			treeFileName = "sim" + nSim + "_" + protocol + ".nwk"

			# Write the tree to file
			tree = Phylo.read(StringIO(tree), 'newick')
			Phylo.write(tree, directory + treeFileName, 'newick')


	##################################################
	## Get migration events counts and rates 
	counts = countsAndRates(t, nDays)

	##################################################
	## Make the last modifications
	if protocol == "all":
		out = counts
		out['nSeq'] = nSeq

	else:
		supp = pd.DataFrame({'value':[rootSDate, rootSLocation, treeHeight], 
			'parameter':['mrcaDate', 'mrcaLocation', 'treeHeight']})

		out = pd.concat([counts, supp], ignore_index=True, sort=False)
		out['nSeq'] = int(nSeq)
		

	## Add columns
	out['protocol'] = re.sub("_[0-9]*$", "", protocol)
	out['matrix'] = matrix
	out['nSim'] = int(nSim)

	print(filename, protocol, 'Done')

	return(out)










def writeNewickTree(filename, transChain, protocol, directory = None):
	"""
	Function that writes the phylogenetic tree in Newick format. 
	It is largely based on the migrationEvents function.

	Arguments:
		- filename (str): name of the simulated transmission chain
		- transChain (dataframe): simulated transmission chain
		- protocol (str): name of the sampling protocol 
			- "all": considers the complete phylogeny and doesn't write a Newick tree
			- see all protocols in the "protocols" variable of the simulatedTreesData module
		- directory (str): name of the directory where to write the Newick tree file  


	Outputs: 
		- It writes the Newick tree when a specific sampling protocol is passed.
		The file name is simX_protocol_nSeq.nwk 
	"""

	##################################################
	## Output directory name 
	if directory:
		if not re.search(r'/$', directory):
			directory += "/"
	else:
		directory = ""

	##################################################
	## Get information from the name of the file
	# nSim 
	nSim = re.findall(r'sim([0-9]+)_', filename)[0]

	##################################################
	## Verify that the sampling protocol is correct
	if protocol not in protocols:
		raise ValueError("The value of the protocol argument is not referenced. \
			\nReferenced values are 'all' and values in the simulatedTreesData module.")

	##################################################
	## Get the list of all the common ancestors of 
	## individuals that were sampled
	# Subset the dataframe
	cols = [
		'case', 
		'source',
		'samp.dates',
		protocol
		]
	transChain = transChain[cols] 

	# List of individuals in the sample
	sample = transChain['case'][transChain[protocol] == 1]
	sample = list(sample)

	## Sampling date of the last sample
	lastDate = transChain['samp.dates'][transChain[protocol] == 1].sort_values(ascending = False).iloc[0]
	
	## Get all the ancestors 
	ancestors = [getAncestors(x, transChain) for x in sample]

	## Intersection of ancestors 
	inter = set.intersection(*map(set, ancestors))
	interL = list(inter)

	## Identification of sampling date and location of MRCA
	sub = transChain[transChain['case'].isin(interL)].sort_values(by='samp.dates', ascending = False)
	rootS = sub['case'].iloc[0]
	rootSDate = sub['samp.dates'].iloc[0]

	# nSeq
	nSeq = re.findall("[150]+", protocol)[0]
	if int(nSeq) not in sampleSize:
		raise ValueError("{0} is not correctly formated.".format(protocol))

	## Trim phylogeny
	t = transChain[(transChain['samp.dates'] >= rootSDate) & (transChain['samp.dates'] <= lastDate)].copy().reset_index()

	##################################################
	## Create the Newick tree
	df = trimPhylogeny2(t, sample, ancestors, inter, rootS)
	tree = newickTree(rootS, df, rootSDate)

	# Newick tree file name
	treeFileName = "sim" + nSim + "_" + protocol + ".nwk"

	# Write the tree to file
	tree = Phylo.read(StringIO(tree), 'newick')
	Phylo.write(tree, directory + treeFileName, 'newick')










def associationIndex(treeFileName, traitFileName, nShuffles, 
	directory = None, seed = None): #, parallelize = False):
	""""""
	print(treeFileName)

	# Directory
	directory = checkDirectory(directory)

	# Read tree_test
	tree = dendropy.Tree.get(file=open(directory + treeFileName, "r"), schema="newick")

	# Read regions file
	regions = pd.read_csv(directory + traitFileName, header = 0, sep = "\t")
	#regions['tips'] = regions.traits.str.replace('_[0-9.]+$', '')
	regions.rename(columns={'traits':'tips'}, inplace=True)

	# Randomly shuffle tip locations
	locations = []

	if seed:
		random.seed(seed)

	for i in range(nShuffles):
		col = 'random' + str(i+1)
		regions[col] = regions.regions.sample(frac = 1).reset_index(drop=True)
		locations.append(regions[['tips', col]])

	# Compute d_AI for the trees with shuffled locations
	random_ai = []
	#if parallelize:

	# def helperF(l):
	# 	return(d_associationIndex(tree, l))

	# with ProcessPoolExecutor() as executor:
	# 	for result in executor.map(helperF, locations):
	# 		random_ai.append(result)

	# else:
	for chars in locations:
		random_ai.append(d_associationIndex(tree, chars))

	# Compute d_AI for the real tree
	tree_ai = d_associationIndex(tree, regions[['tips', 'regions']])

	# Compute the association index
	AI = tree_ai / mean(random_ai)
	return(AI)










def d_associationIndex(tree, chars):
	"""

	"""
	# Rename the column corresponding to the locations
	chars = chars.rename(columns=lambda s: 'regions' if not s.find('random') else s)

	d_ai = []

	for node in tree.preorder_internal_node_iter():

		# Get the tips of the subtree whose root is node 
		if node.parent_node is None:

			# Number of traits samples
			c = len(chars.regions.unique())

			# Frequency table of each trait
			freq_tab = chars.groupby('regions').size() / chars.shape[0]

		else:

			# List of tips
			tips_nodes = node.leaf_nodes()
			tips = [re.sub("'", "", str(x.taxon)) for x in tips_nodes]
			tips = [x.replace(" ", "_") for x in tips]
			
			# Number of traits samples
			subChars = chars[chars.tips.isin(tips)]
			c = len(subChars.regions.unique())

			# Frequency table
			freq_tab = subChars.groupby('regions').size() / subChars.shape[0]


		f = max(freq_tab)
		s = (1 - f) / (2*c - 1)
		d_ai.append(s)

	return(sum(d_ai))










def migrationEffectiveness(fileName, directory, matrix):
	"""

	"""
	
	##################################################
	## Directory name 
	if directory:
		if not re.search(r'/$', directory):
			directory += "/"
	else:
		directory = ""

	##################################################
	## Load and rearrange the transmission chain 
	# Load
	transChain = pd.read_csv(directory + fileName, sep="\t", header = 0)

	# Add regions.source column to transChain
	transChain['patch'] = transChain.case.str.replace(r'_.*$', '').astype('int')
	transChain['patch'] += 1

	transChain['patch.source'] = transChain.source.str.replace(r'_.*$', '').astype('int')
	transChain['patch.source'] += 1

	transChain['regions.source'] = transChain['patch.source'].replace(regionDic)
	transChain['regions'] = transChain['patch'].replace(regionDic)

	# Subset dataframe 
	transChain = transChain[['case', 'source', 'regions', 'regions.source']]

	# nSim 
	nSim = re.findall(r'sim([0-9]+)_', fileName)[0]

	##################################################
	## Identify migration events
	migrants = transChain[transChain['regions'] != transChain['regions.source']]

	## Apply getTransmissionChain function
	out = migrants.apply(transmissionChainLength, axis = 1, transChain = transChain)

	out.reset_index(inplace = True, drop = True)
		
	##################################################
	## Add columns
	out['matrix'] = matrix
	out['nSim'] = int(nSim)

	print(fileName, 'Done')

	return(out)










def transmissionChainLength(migrant, transChain):
	# Retrieve arguments 
	indexCase = migrant['source']
	case = migrant['case']
	source = migrant['regions.source']
	destination = migrant['regions']

	# Trim the transmission chain to the destination region
	df = transChain.copy()
	df = df[(df['regions'] == destination) & (df['regions.source'] == destination)]

	# Get the transmission list
	transChainList = getTransmissionChain(case, df)
	length = len(transChainList) + 1

	out = pd.Series({'source': source, 
		'destination': destination, 
		'length' : length})

	return(out)
	









def getTransmissionChain(case, df):
	""""""
	# Conditions and outputs
	descCond = (df['source'] == case) # Condition to get the offspring of ind
	descN = sum(descCond)
	desc = list(df[descCond]['case']) # List of case's offspring

	# ind has at least one descendant
	if descN != 0:
		if descN == 1: # One descendant
			out = []
			out.append(desc)
			innerDesc = getTransmissionChain(desc[0], df)
			if len(innerDesc) == 1:
				out.append(innerDesc)
			else:
				out += innerDesc
			return out

		else: # Several descendants
			out = []
			out.append(desc)
			for d in desc:
				innerDesc = getTransmissionChain(d, df)
				if len(innerDesc) == 1 :
					out.append(innerDesc)
				else:
					out += innerDesc
				return out

	# ind is a tip of the tree
	else:
		out = []
		return(out)





