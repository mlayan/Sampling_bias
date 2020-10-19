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
import dendropy
import itertools
import math
import numpy as np
import os
import pandas as pd
import random
import re
import sys

# Import functions
from Bio import Phylo
from concurrent.futures import ProcessPoolExecutor
from io import StringIO
from itertools import permutations
from statistics import mean

# Import fixed data from trees
from handleDirectories import checkDirectory
from simulatedTreesData import *








 
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

	while not ancestor in indexCase:
		ancestor = transmissionChain.loc[transmissionChain['case'] == individual, 'source'].values[0]
		ancestorsList.append(ancestor)
		individual = ancestor

	return(ancestorsList)










def countsAndRates(df, t, regionDic):
	""" 
	Function to compute the number of migration events from a transmission chain dataframe
	Parameters:
		- df (pandas dataframe): transmission chain with columns
			regions : names of cases' region
			regions.source : names of sources' region
		- time (int) : time span of the transmission chain in days
		- regionDic (dict): dictionary of the identifiers (keys) and the region names (values)

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

	# Migration list
	migrationList = list(regionDic.values())
	migrationList = list(itertools.permutations(migrationList, 2))
	migrationList = [x + "_" + y for x,y in migrationList]

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










def reconstructSampleName(x):
	"""
	"""
	if isinstance(x, np.integer):
		return(str(math.ceil(x)))
	elif x.is_integer(): 
		return(str(math.ceil(x)))
	else: 
		return(str(round(x, 6)))










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
	
	# Rename columns
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

		labelString = reconstructSampleName(label[0])

		out = ind + "_" + labelString
		return(out, label[0])










def newickTree(rootTree, df, rootDate, innerlabel = True):
	"""
	Function that returns the final Newick tree from a phylogenetic tree in a dataframe format.
	It is based on the treeScan(NoLabel) functions.
	Arguments:
		- rootTree (string): name of the tree root
		- df (dataframe) : phylogenetic tree
			- column 1: nodes name
			- column 2: children name
			- column 3: sampling date of children
		- rootDate (int): sampling date of rootTree
		- innerlabel (bool): add branch length in the Newick tree

	Output:
		- out (string): Final Newick tree
	"""
	if innerlabel:
		(out, nodeDate) = treeScan(ind = rootTree, df = df, nodeDate = rootDate) 

		if nodeDate != rootDate:
			out = "(" + out 
			out += ":" + str(nodeDate - rootDate) + ");"
		else:
			out += ";"

		return(out)

	else:
		out = treeScanNoLabel(rootTree, df)
		out += ";"
		return(out)










def nexusTree(rootTree, df, rootDate, rootLocation, sampleCompleteNames, filename):
	"""
	Function that returns the final Newick tree from a phylogenetic tree in a dataframe format.
	It is based on the treeScan(NoLabel) functions.
	Arguments:
		- rootTree (string): name of the tree root
		- df (dataframe) : phylogenetic tree
			- column 1: nodes name
			- column 2: children name
			- column 3: sampling date of children
			- column 4: children location
			- column 5: nodes location
		- rootDate (int): sampling date of rootTree
		- innerlabel (bool): add branch length in the Newick tree

	Output:
		- out (string): Final Newick tree
	"""
	if len(df.columns) != 5:
		raise ValueError("Dataframe has not the correct dimension : {} columns instead of 5".format(len(df.columns)))

	## Change column names
	df.columns = ['source', 'child', 'samp.dates', 'regions', 'regions.source']

	# Assign rootDate when not specified
	if not rootDate : 
		rootDate = df[df['child'] == ind]['label'].values
		if len(rootDate) > 1:
			raise ValueError("{0} appears multiple times in the dataframe as a tip".format(ind))
		rootDate = rootDate[0]

	# Create empty tree with taxon names 
	taxon_namespace = dendropy.TaxonNamespace(sampleCompleteNames)
	tree = dendropy.Tree(taxon_namespace = taxon_namespace)
 
	# Add nodes and annotations programmatically
	out = treeScanNexus(tree, df, rootTree, sourceDate = rootDate, 
		sourceLoc = rootLocation, treeRoot = True, filename = filename)

	return(out)









def treeScanNexus(tree, df, ind, nMig = None, sourceDate = None, sourceLoc = None, 
                  treeRoot = False, filename = None, **detailedCounts):
	"""
	Scan the simulated transmission chain to write a tree with complete annotations
	"""
	# Offspring of ind
	descCond = (df['source'] == ind)    # Condition to get the offspring of ind
	desc = list(df[descCond]['child'])  # List of ind's offspring
	descN = sum(descCond)               # Size of ind's offspring

	# ind has at least one descendant
	if descN != 0:

		# Ind characteristics
		if treeRoot: # For the tree root
			indLoc = sourceLoc
			indDate = sourceDate
			# Create seed node
			tree.seed_node.regions_states = indLoc
			tree.seed_node.annotations.add_bound_attribute("regions_states", 
				annotation_name="regions.states")

		else:        # For inner nodes
			indLoc = df[descCond]['regions.source'].iloc[0]
			indDate = df[df['child'] == ind]['samp.dates'].iloc[0]
			sourceLoc = df[df['child'] == ind]['regions.source'].iloc[0]
			# Migration events
			if not nMig: nMig = 0
			if indLoc != sourceLoc: nMig += 1
			# Detailed migration events
			if detailedCounts:
				detailedCountsInner = getDetailedCounts(indLoc, sourceLoc, df, **detailedCounts)
			else:
				detailedCountsInner = getDetailedCounts(indLoc, sourceLoc, df)


		if descN == 1: # One descendant
			if treeRoot:
				childNode = treeScanNexus(tree, df, desc[0], sourceDate=indDate)
				tree.seed_node.add_child(childNode)
				return(tree)
			else:
				return treeScanNexus(tree, df, desc[0], nMig = nMig, 
					sourceDate = sourceDate, **detailedCountsInner)               

		else: # Several descendants

			if treeRoot:
				# Add children nodes
				for d in desc:
					childNode = treeScanNexus(tree, df, d, sourceDate = indDate)
					tree.seed_node.add_child(childNode)
				return(tree)

			else:
				# Create inner node
				innerNode = dendropy.Node(edge_length = indDate - sourceDate)
				innerNode = addDetailedCountsToNode(innerNode, indLoc, sourceLoc, 
					df, nMig, **detailedCountsInner)

				# Add children nodes
				for d in desc:
					childNode = treeScanNexus(tree, df, d, sourceDate = indDate)
					innerNode.add_child(childNode)

				return(innerNode)


	# ind is a tip of the tree
	else:
		# Ind characteristics
		indDate = df[df['child'] == ind]['samp.dates'].values[0]
		indLoc = df[df['child'] == ind]['regions'].values[0]
		sourceLoc = df[df['child'] == ind]['regions.source'].values[0]

		# Number of migrations
		if not nMig: nMig = 0
		if indLoc != sourceLoc: nMig += 1
		if detailedCounts:
			detailedCountsTip = getDetailedCounts(indLoc, sourceLoc, df, **detailedCounts)
		else:
			detailedCountsTip = getDetailedCounts(indLoc, sourceLoc, df)   

		# Create dendropy node 
		tipNode = dendropy.Node(edge_length = indDate - sourceDate)
		tipNode = addDetailedCountsToNode(tipNode, indLoc, sourceLoc, df, nMig, **detailedCountsTip)

		# Taxon name 
		indDateString = reconstructSampleName(indDate)
		tipNode.taxon = tree.taxon_namespace.get_taxon(ind + "_" + indDateString)

		return(tipNode)











def addDetailedCountsToNode(node, indLoc, sourceLoc, df, nMig, **detailedCounts):
	"""
	Add migration counts to the substanding edge of the node  
	"""    
	# # Get the dictionary of the detailed counts
	# regionCounts = getDetailedCounts(indLoc, sourceLoc, df, **detailedCounts)
	# Add the detailed counts to the node
	for k, v in detailedCounts.items():
		node.edge.annotations.add_new(k, v)
		#node.edge.annotations.add_bound_attribute(k) 

	# Add total counts
	node.edge.regions_count = nMig
	node.edge.annotations.add_bound_attribute("regions_count", 
		annotation_name="regions.count")

	# Add region states
	node.regions_states = indLoc
	node.annotations.add_bound_attribute("regions_states", 
		annotation_name="regions.states")

	return(node)










def getDetailedCounts(indLoc, sourceLoc, df, **detailedCounts):
	"""
	Create or edit a dictionary of migration counts 
	"""
	# Create a dictionary if it doesn't exist 
	if not detailedCounts:
		# Get the set of regions if detailedCounts does not exist
		regionSet = set(list(df['regions'].unique()) + list(df['regions.source'].unique()))
		regionSet = permutations(regionSet, 2)
		regionPermutations = ["counts_" + x + "_" + y for x,y in regionSet]
		# Create "empty" dictionary
		detailedCounts = {k:0 for k in regionPermutations}

	# Update if necessary the dictionary
	if indLoc != sourceLoc and sourceLoc:
		detailedCounts['counts_' + sourceLoc + "_" + indLoc] += 1

	return(detailedCounts)










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
			nodeDate = out[out.case == x ]['samp.dates'].iloc[0]        
			out.replace({x:addNode}, inplace = True)
			out = out.append({'case': x, 'source': addNode, 'samp.dates': nodeDate + 0.001}, 
				ignore_index = True)

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
	out = transChain[['source', 'case', 'samp.dates', 'regions', 'regions.source']]
	out = out[out.source.isin(ancestorsSubtree)]
	out = out[out.case.isin(ancestorsSubtree + tips)]

	# Additionnal nodes for sample individuals which are also part 
	## of the inner phylogeny
	for x in tips:
		if x in ancestorsSubtree:
			addNode = "a_" + x    # Additionnal node which will be considered as the ancestor of x
			nodeDate = out[out['case'] == x ]['samp.dates'].iloc[0]   # Sampling date of x
			nodeLocation = out[out.case == x]['regions'].iloc[0] 	# Sampling location of x
			nodeSourceLocation = out[out.case == x]['regions.source'].iloc[0] 	# Sampling location of the source of x
			out.at[out['case'] == x, 'samp.dates'] = nodeDate - 0.000001
			out.replace({x:addNode}, inplace = True)   # All occurences of x are replaced by the additionnal node
			out = out.append({
				'case':x, 
				'source':addNode, 
				'samp.dates':nodeDate, 
				'regions':nodeLocation, 
				'regions.source':nodeSourceLocation}, 
			ignore_index = True) # new row corresponding to the branch (additionnalnode -> x)

	return(out)










def migrationEvents(filename, matrix, protocol, regionDic, 
	extractNewickTree = False, extractNexusTree = False, directory = None):
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
	if extractNewickTree and protocol == "all":
		raise ValueError('No tree is extracted when migration events \
			are evaluated over the entire phylogeny.')

	##################################################
	## Directory name 
	directory = checkDirectory(directory)

	##################################################
	## Load and rearrange the transmission chain 
	# Load
	transChain = pd.read_csv(directory + filename, sep="\t", header = 0)

	# Drop the row corresponding to the root
	transChain = transChain[transChain['case'] != indexCase[0]]

	# Add regions.source column to transChain
	transChain['regions.source'] = transChain['patch.source'].replace(regionDic)

	# nSim 
	nSim = re.findall(r'sim([0-9]+)_', filename)[0]

	# Location of the start of the epidemic
	startLocation = list(regionDic.values())[0]
	
	##################################################
	## Extract informations 
	if protocol == "all":
		# Number of days corresponding to the entire epidemic
		nDays = round(365.25*nYears)

		# nSeq doesn't take the size of the sample but na
		nSeq = float("nan")

	else:
		if protocol not in protocols:
			raise ValueError("The value of the protocol argument is not referenced. \
				\nReferenced values are 'all' and values in the simulatedTreesData module.")

		##################################################
		## Get the list of all the common ancestors of 
		## individuals that were sampled
		# List of individuals in the sample
		sample = transChain['case'][transChain[protocol] == 1]
		sampleCompleteNames = sample + "_" + \
							transChain[transChain[protocol] == 1]['samp.dates'].apply(lambda x: reconstructSampleName(x))   
		sample = list(sample)

		## Sampling date of the last sample
		lastDate = transChain['samp.dates'][transChain[protocol] == 1].sort_values(ascending = False).iloc[0]
		
		## Get all the ancestors 
		ancestors = [getAncestors(x, transChain[['case', 'source']]) for x in sample]

		## Intersection of ancestors 
		inter = set.intersection(*map(set, ancestors))
		interL = list(inter)
		if len(interL) == 1:
			if indexCase == interL:
				root = indexCase[0]
				rootLocation = startLocation
				rootDate = startDate
			else:
				raise ValueError("A single MCRA has been found : {0}. \
					It does not correspond to {1}.".format(interL, indexCase))
		else:
			## Identification of sampling date and location of MRCA
			sub = transChain[transChain['case'].isin(interL)].sort_values(by='samp.dates', ascending = False)
			root = sub['case'].iloc[0]
			rootLocation = sub['regions'].iloc[0]
			rootDate = sub['samp.dates'].iloc[0]

		# nSeq
		nSeq = re.match(r'.*_(\d+)$', protocol).group(1)
		if int(nSeq) not in sampleSize:
			raise ValueError("{0} is not correctly formated.".format(protocol))

		# EpidemicDuration in days  
		treeHeight = lastDate - rootDate
		nDays = round(treeHeight*365.35)

		## Trim phylogeny
		t = trimPhylogeny2(transChain, sample, ancestors, inter, root)
		# t = transChain[(transChain['samp.dates'] >= rootDate) & (transChain['samp.dates'] <= lastDate)].copy().reset_index()

		##################################################
		# Write nexus tree
		if extractNexusTree:
			tree = nexusTree(root, t, rootDate, rootLocation, sampleCompleteNames, filename)
			treeFileName = "sim" + nSim + "_" + protocol + ".nex"
			tree.write(path=directory + treeFileName, schema='nexus')

		# Write Newick tree
		if extractNewickTree:
			tree = newickTree(root, t, rootDate)
			treeFileName = "sim" + nSim + "_" + protocol + ".nwk"
			tree = Phylo.read(StringIO(tree), 'newick')
			Phylo.write(tree, directory + treeFileName, 'newick')

	##################################################
	## Get migration events counts and rates 
	if protocol == "all":
		counts = countsAndRates(transChain, nDays, regionDic)
	else:
		counts = countsAndRates(t, nDays, regionDic)

	##################################################
	## Make the last modifications
	if protocol == "all":
		out = counts
		out['nSeq'] = nSeq

	else:
		supp = pd.DataFrame({'value':[rootDate, rootLocation, treeHeight], 
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
	root = sub['case'].iloc[0]
	rootDate = sub['samp.dates'].iloc[0]

	# nSeq
	nSeq = re.findall("[150]+", protocol)[0]
	if int(nSeq) not in sampleSize:
		raise ValueError("{0} is not correctly formated.".format(protocol))

	## Trim phylogeny
	t = transChain[(transChain['samp.dates'] >= rootDate) & (transChain['samp.dates'] <= lastDate)].copy().reset_index()

	##################################################
	## Create the Newick tree
	df = trimPhylogeny2(t, sample, ancestors, inter, root)
	tree = newickTree(root, df, rootDate)

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










def migrationEffectiveness(fileName, directory, matrix, 
	regionDic):
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





