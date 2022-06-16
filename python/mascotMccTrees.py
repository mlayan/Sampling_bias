#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MODULE TO COMPUTE MASCOT MCC TREES FROM MAPPED
TREE POSTERIOR DISTRIBUTION
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2022-01-17' 
__last_update__ = ''

# Import libraries
import os
import sys
import re
import math
import dendropy as dp
import pandas as pd

## Import custom functions
from handleDirectories import checkDirectory







def removeUnifurcations(tree):
	"""
	Remove all unifurcations (nodes with outdegree 1) from dendropy.Tree object
	"""
	filter_fn = lambda n: (len(n.child_nodes())!=1 or n.parent_node is None) and not n.is_leaf() 

	for nd in tree.preorder_node_iter(filter_fn = filter_fn):   
		for child in nd.child_nodes():
			grandchildren = child.child_nodes() 
			while len(grandchildren) == 1:
				grandchildren[0].edge_length += child.edge_length
				grandchildren[0].parent_node = nd
				nd.remove_child(child, suppress_unifurcations=False)
				nd.add_child(grandchildren[0])
				child = grandchildren[0]
				grandchildren = child.child_nodes()

	out = tree.clone()
	out.purge_taxon_namespace()
	out.encode_bipartitions(suppress_unifurcations=False, collapse_unrooted_basal_bifurcation=False)

	return(out)










def mccTree(fileName, directory, burnIn=10):
	"""
	Get MCC tree for a tree posterior distribution obtained by MASCOT or MASCOT-GLM 
	"""

	print(fileName)
	
	# Burn-in
	directory = checkDirectory(directory)
	logFile = pd.read_csv(directory+fileName.replace('trees', "log"), comment='#', sep="\t")
	treeFiles = [directory + fileName]

	if os.path.exists(directory+fileName.replace('.trees', "_b.log")):
		temp_file = pd.read_csv(directory+fileName.replace('.trees', "_b.log"), comment='#', sep="\t")
		logFile = pd.concat([logFile, temp_file], sort = False, ignore_index = True)
		treeFiles.append(directory + fileName.replace('.trees', '_b.trees'))
	b = math.ceil(logFile.shape[0] * burnIn / 100)

	# Iterate over post burn-in posterior distribution and store trees	
	if "mascot" in directory or "glm" in directory:
		extract_node_info = True
	else :
		extract_node_info = False

	treeCollectionYielder = dp.Tree.yield_from_files(
		files = treeFiles, 
		schema = 'nexus',  
		preserve_underscores = True,  
		extract_comment_metadata = extract_node_info,
		rooting = 'default-rooted'
	)

	treeList = dp.TreeList()
	for i, tree in enumerate(treeCollectionYielder):
		if i > b:
			if "mascot" in directory or "glm" in directory:
				treeList.append(removeUnifurcations(tree))
			else:
				treeList.append(tree)

	# Get and write mcc tree
	mcc = treeList.maximum_product_of_split_support_tree()
	mcc.write_to_path(directory + fileName.replace("trees.txt", "mcc.tree"), "nexus")

	return(0)

