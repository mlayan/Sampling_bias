#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MODULE TO COMPARE TREE TOPOLOGIES
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-05-25' 
__last_update__ = '2020-05-'

# Import libraries
import os
import sys
import re
import pandas as pd
import numpy as np
import dendropy

# Import functions
from itertools import combinations
from sklearn.linear_model import LinearRegression

# Custom functions
from handleDirectories import checkDirectory










def comparePairwisetMRCA(fileName, sim, **mcc):
	"""
	Function to compare pairwaise times to the MRCA between the mcc tree 
	and the simulated tree 
		- sim (dendropy tree) : simulated tree
		- **mcc (dendropy trees) : estimated mcc trees
		keys should corrrespond to the estimation model

	It returns a dataframe with the multiple informations:
		- inference model 
		- R2 
		- intercept
		- slope
	"""

	# Verify that lists of taxa are identical
	simTaxa = sorted([x.label for x in sim.taxon_namespace])

	for k, v in mcc.items():
		mccTaxa = sorted([x.label for x in v.taxon_namespace])
		if mccTaxa != simTaxa:
			# print('simulation' + str(nSim) + '/' + k + '/' + fileName)
			raise ValueError('{} and simulated trees do not have the same taxa. See file {}'\
	 		 	.format(k, fileName))

	# Get patristic distances into a list of tuples
	trees = mcc.copy()
	trees['sim'] = sim

	for i, (k, v) in enumerate(trees.items()):
		
		# Patristic distance matrix 
		dist = v.phylogenetic_distance_matrix().as_data_table()

		# Store in a list of tuples
		distances = []
		for row in dist.row_name_iter():
			for col in dist.column_name_iter():  
				if col != row:
					distances.append((row, col, dist[row, col]))
		
		# Coerce into a dataframe
		distData = pd.DataFrame(distances, columns = ['row', 'col', k])

		if i == 0:
			merged = distData.copy()
		else:
			merged = pd.merge(merged, distData, on = ['row', 'col'])


	# Convert dataframe columns into numpy arrays and 
	# perform linear regression
	coef = []
	intercept = []
	slope = []

	X = merged.sim.to_numpy().reshape(-1,1)

	for m in mcc.keys():
		# Convert to numpy array
		Y = merged[m].to_numpy().reshape(-1, 1)

		# Linear regression
		linear_regressor = LinearRegression()
		lm = linear_regressor.fit(X, Y)

		# Parameters 
		slope += lm.coef_[0].tolist()
		intercept += lm.intercept_.tolist()
		coef.append(lm.score(X, Y))

	out = pd.DataFrame({'model' : list(mcc.keys()), 
		'coef': coef,
		'intercept' : intercept,
		'slope' : slope})

	return(out)











def sort_namespace(namespace):
	"""
	Function which sorts a taxon namespace object 
	(dendropy.datamodel.taxonmodel.TaxonNamespace) according to its 
	labels. 

	The function doesn't return a taxon namespace object but
	a list of taxon objects ordered according to their label.
	"""
	labs = sorted(namespace.labels())
	labsDict = namespace.label_taxon_map()

	out = []
	for lab in labs:
		out.append(labsDict.get(lab))

	return(out)










def compareBeasttoSimulation(mccFile, matrix, mccDirectory, *models):
	"""
	Wrapper function which computes the linear regression of pairwise TMCRA between 
	simulated and inferred phylogenetic trees 
		- mccFile (str): name of the mcc tree file
		- matrix (int) : number of the spatial coupling matrix
		- mccDirectory (str): directory where to find the mcc tree file
		- models (list of strings): list of Beast models for which 
		the tree topology is compared to the real tree topology

	It returns a pandas dataframe with the following columns:
		- number of simulation
		- name of the sampling protocol
		- sample size
		- number of the spatial coupling matrix
		- Beast inference model
		- regression coefficient
		- slope of the linear regression 
		- intercept of the linear regression  
	"""
	
	print(mccFile)

	# Get general information
	nSim = re.match(r'sim(\d+).*', mccFile).group(1)
	nSeq = re.match(r'.*_(\d+)\.mcc.*', mccFile).group(1)
	protocol = re.match(r'sim\d+_(.*)_.*', mccFile).group(1)
	
	# Correct directory
	mccDirectory = checkDirectory(mccDirectory)
	directories = {}
	to_replace = re.match(r'.*/(\w+)/$', mccDirectory).group(1)

	for m in models:
		if m in mccDirectory:
			directories[m] = mccDirectory
		else:
			directories[m] = mccDirectory.replace(to_replace, m)

	# Get simulated tree
	simDirectory = mccDirectory.replace(to_replace, 'files')
	simFile = mccFile.replace('mcc.tree', 'nex')
	sim = dendropy.Tree.get(path = simDirectory + simFile, 
			schema = 'nexus', 
			preserve_underscores = True)

	# Get inferred mcc trees
	trees = {}
	for k, v in directories.items():

		if os.path.exists(v + mccFile):
			trees[k] = dendropy.Tree.get(path = v + mccFile, 
				schema = 'nexus', 
				preserve_underscores = True)

	# Get linear regression coefficient
	data = comparePairwisetMRCA(mccFile, sim, **trees)

	# Additionnal information
	data['nSim'] = nSim
	data['protocol'] = protocol
	data['nSeq'] = nSeq
	data['matrix']  = matrix

	return(data)