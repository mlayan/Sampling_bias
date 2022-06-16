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

## Import custom functions
from handleDirectories import checkDirectory
from mathFunctions import *
from mascotMccTrees import removeUnifurcations

## Import regions and log files columns dictionaries
from beastLogDictionaries import * 

# Do not iterate over nodes with a single child (applicable to mascot)
# except if it is the root 
filter_fn = lambda n: len(n.child_nodes()) != 1 or n.parent_node == None  



def treePosteriorDistributionStatistics(fileName, directory, runType, regions, 
	lineage_migration = True, introduction_dates = True, burnIn=10):

	print(fileName)
	post_dist = treePosteriorDistribution(fileName, directory, runType, regions, burnIn)
	post_dist.calc_statistics(lineage_migration, introduction_dates)
	lineage_migeation_out, introduction_dates_out = post_dist.get_statistics()

	return lineage_migeation_out, introduction_dates_out



def summary_statistics(posterior_numpy_array, statistic, regions):
	
	# Percentiles for WIS calculation
	all_percentiles= [x for x in range(5, 100, 5)] + [1, 2.5, 99, 97.5]
	all_percentiles = sorted([x/100 for x in all_percentiles])

	# Compute summary statistics
	all_statistics = []

	all_statistics.append(np.apply_along_axis(mean, 0, posterior_numpy_array)) # Mean

	all_statistics.append(np.apply_along_axis(									# lower bound of 95% HPD
		lambda x: hpd(x, lower_only = True), 0, posterior_numpy_array)
	)

	all_statistics.append(np.apply_along_axis(									# upper bound of 95% HPD
		lambda x: hpd(x, upper_only = True), 0, posterior_numpy_array)
	)

	for q in all_percentiles:
		all_statistics.append(np.quantile(posterior_numpy_array, q, axis = 0))		# quantiles

	# Store in a pandas dataframe
	statistics = ['mean', 'X2_5_hpd', 'X97_5_hpd']
	statistics += ['X'+str(round(q*100,1)).replace(".0", "").replace(".5", "_5") for q in all_percentiles]

	if statistic == "introduction_dates":
		out = pd.DataFrame({"regions": list(regions.keys())})
		for i, m in enumerate(all_statistics):
			tempDF = pd.DataFrame({statistics[i] : list(all_statistics[i])}) 
			out = pd.concat([out, tempDF], axis = 1)

	if statistic == "lineage_migration":
		out = pd.DataFrame()
		for i, m in enumerate(all_statistics):
			tempDF = pd.DataFrame(m, index = list(regions.keys()), columns = list(regions.keys())).reset_index().rename(
				columns = {'index':'source'}
			)
			tempDF = pd.melt(tempDF, id_vars=['source'], value_vars=list(regions.keys()), 
				var_name='destination', value_name=statistics[i])    
			if out.empty:
				out = tempDF
			else:
				out = pd.merge(out, tempDF, on = ['source', 'destination'], sort = True)

	return out	



class treePosteriorDistribution:

	def __init__(self, fileName, directory, runType, regions, burnIn):

		# Annotation name
		beastModel = runType
		if 'glm' in runType: beastModel = 'mascot'
		if 'files' in runType: beastModel = 'simulation'
		self.ann = annotationDict[beastModel]

		# Basic attributes
		directory = checkDirectory(directory)
		self.directory = directory
		self.name = fileName
		self.regions = regions
		self.nregions = len(regions)
		self.model = runType
		self.nSim = int(re.match(r'^.*sim(\d*)_', fileName).group(1))
		self.nSeq = int(re.match(r'^.*_(\d*)\.[nextrees]+', fileName).group(1))
		self.protocol = re.match(r'^.*sim\d*_(.*)_\d*\.[nextrees]+', fileName).group(1)

		self.lineage_migration = None
		self.introduction_dates = None

		# Burn-in
		if runType == "files":
			self.burnin = None
			self.treeFiles = None
			self.simtree = directory+fileName

		else :
			logName = fileName.replace("trees", "log")
			multipleBastaChain = len([f for f in os.listdir(directory) if re.match(fileName.replace(".trees.txt", "")+".*trees.txt", f)]) > 1 and beastModel == "basta"
			resumedMascotChain = os.path.exists(directory+fileName.replace(".trees", "_b.trees")) and beastModel == "mascot"

			if not burnIn and not multipleBastaChain: 
				burnIn = [0]
			elif not burnIn and multipleBastaChain:
				burnIn = [0] * len([f for f in os.listdir(directory) if re.match(fileName.replace(".trees.txt", "")+".*trees.txt", f)])
			elif resumedMascotChain:
				logFile = pd.concat(
					[
						pd.read_csv(directory + logName, sep = "\t", comment = "#"), 
						pd.read_csv(directory + logName.replace('.log', '_b.log'), sep = "\t", comment = "#")
					],
					sort = False, ignore_index = True) 
				chainLength = logFile.shape[0] 
				burnIn = [round(chainLength * burnIn / 100)]
			elif multipleBastaChain:
				b = burnIn
				burnIn = []
				for f in [directory+x.replace("trees", "log") for x in os.listdir(directory) if re.match(fileName.replace(".trees.txt", "")+".*trees.txt", x)]:
					logFile = pd.read_csv(f, sep = "\t", comment = "#")
					chainLength = logFile.shape[0] 
					burnIn.append(round(chainLength * b / 100))  
			else:
				logFile = pd.read_csv(directory + logName, sep = "\t", comment = "#")
				chainLength = logFile.shape[0] 
				burnIn = [round(chainLength * burnIn / 100)]

			# Create tree yielder to avoid memory outspacing
			treeFiles = []
			if multipleBastaChain:
				treeFiles = [[directory+f] for f in os.listdir(directory) if re.match(fileName.replace(".trees.txt", "")+".*trees.txt", f)]
			elif resumedMascotChain:
				treeFiles = [[directory+fileName, directory + fileName.replace(".trees", "_b.trees")]]
			else :
				treeFiles = [[directory+fileName]]
			

			self.burnin = burnIn
			self.treeFiles = treeFiles
			self.simtree = None 






	def calc_statistics(self, lineage_migration, introduction_dates):

		if self.model == "files":
			return self.calc_statistics_sim(lineage_migration, introduction_dates)

		else:
			return self.calc_statistics_beast(lineage_migration, introduction_dates)






	def calc_statistics_sim(self, lineage_migration, introduction_dates):

		# Simulated transmission chain in nexus format
		tree = dp.Tree.get(path=self.simtree, schema="nexus")

		# Initialize transition matrix and introduction dates
		transitionMatrix = np.zeros((self.nregions, self.nregions))
		
		dateVector = np.zeros(self.nregions)
		seedLocation = tree.seed_node.annotations.get_value(self.ann)
		dateVector[self.regions[seedLocation]] = tree.max_distance_from_root()
		
		# Compute statistics
		for node in tree.preorder_node_iter():
			
			source = node.annotations.get_value(self.ann)
			
			for child in node.child_nodes():
				if child is None:
					print("Child node is Nonetype in file:" + self.name)
				
				destination = child.annotations.get_value(self.ann)
				# Add date of introduction on destination 
				if source != destination:
					if dateVector[self.regions[destination]] == 0.0:
						dateVector[self.regions[destination]] = child.root_distance

				transitionMatrix[self.regions[source], self.regions[destination]] += 1
				
				# Stop iterating over nodes if virus is introduced in all locations
				if not np.any( (dateVector == 0.0 ) ) and lineage_migration == False and introduction_dates == True:
					break


		# Store statistics into pandas dataframes
		if lineage_migration == True: 
			lineage_migration_out = pd.DataFrame(transitionMatrix, index=list(self.regions.keys()), columns=list(self.regions.keys)).reset_index().rename(
				columns={'index':'source'}
				)
			self.lineage_migration = pd.melt(lineage_migration_out, 
				id_vars=['source'], value_vars=list(self.regions.keys()), var_name='destination', value_name='value')

		if introduction_dates:
			# Compute calendar dates
			mostRecentSample = max([float(re.sub(r'^.*_', '', l.label)) for l in tree.taxon_namespace])
			dateVector[dateVector == 0] = np.nan
			dateVector *= -1
			dateVector += mostRecentSample

			# Conversion into a pandas dataframe
			self.introduction_dates = pd.DataFrame({
				'regions': list(self.regions.keys()), 
				'value': list(dateVector)
				})








	def calc_statistics_beast(self, lineage_migration, introduction_dates):

		transitionMatrices = []
		dateArray = []

		for j, treeList in enumerate(self.treeFiles):

			treeCollectionYielder = dp.Tree.yield_from_files(
				files = treeList, 
				schema = 'nexus',  
				preserve_underscores = True,  
				extract_comment_metadata = True
				)

			transitions = np.zeros((1, self.nregions, self.nregions))
			dateVector = np.zeros((1, self.nregions))

			for i, tree in enumerate(treeCollectionYielder):
				if i < self.burnin[j]:
					continue

				# Remove  
				if self.model in ["glm", "mascot"]:
					tree = removeUnifurcations(tree)
				
				# Location and temporal info on tree
				seedLocation = tree.seed_node.annotations.get_value(self.ann)
				rootDistances = tree.calc_node_root_distances()
				
				# Object to modify
				if i - self.burnin[j] == 1:
					dateVectorToModify = dateVector
					transitionsToModify = transitions
				else:
					dateVectorToModify = np.zeros((1,self.nregions)) 
					transitionsToModify = np.zeros((1, self.nregions, self.nregions))

				# Edit root introduction date
				dateVectorToModify[0, self.regions[seedLocation]] = tree.max_distance_from_root()     

				# Scan tree and collect migration info
				for node in tree.preorder_node_iter(filter_fn = filter_fn):
					source = node.annotations.get_value(self.ann)

					for child in node.child_nodes():
						destination = child.annotations.get_value(self.ann)

						if lineage_migration == True:
							transitionsToModify[0, self.regions[source], self.regions[destination]] += 1

						if source != destination and introduction_dates == True:
							if dateVectorToModify[0, self.regions[destination]] == 0.0:
								dateVectorToModify[0, self.regions[destination]] = child.root_distance

						# Stop iterating over nodes if virus is introduced in all locations
						if not np.any( (dateVectorToModify == 0.0 ) ) and introduction_dates == True and lineage_migration == False:
							break

				
				# Compute calendar dates
				mostRecentSample = max([float(re.sub(r'^.*_', '', l.label)) for l in tree.taxon_namespace])
				dateVectorToModify[dateVectorToModify == 0] = np.nan
				dateVectorToModify *= -1
				dateVectorToModify += mostRecentSample

				if i - self.burnin[j] > 1:
					dateVector = np.concatenate((dateVector, dateVectorToModify))
					transitions = np.concatenate((transitions, transitionsToModify))

			# Store posterior distribution of independent chains in one vector
			if introduction_dates == True: dateArray.append(dateVector)
			if lineage_migration == True: transitionMatrices.append(transitions)


		# Concatenate list of numpy arrays into a multimensional numpy array
		transitionMatrices = np.concatenate(transitionMatrices, axis = 0)
		dateArray = np.concatenate(dateArray, axis = 0)

		# Store info in treePosteriorDistribution object
		if lineage_migration == True:
			self.lineage_migration = summary_statistics(transitionMatrices, "lineage_migration", self.regions)
		if introduction_dates == True:
			self.introduction_dates = summary_statistics(dateArray, "introduction_dates", self.regions)











	def get_statistics(self):

		if self.lineage_migration is not None:
			self.lineage_migration['model'] = self.model
			self.lineage_migration['nSim'] = self.nSim
			self.lineage_migration['nSeq'] = self.nSeq
			self.lineage_migration['protocol'] = self.protocol

		if self.introduction_dates is not None:
			self.introduction_dates['model'] = self.model
			self.introduction_dates['nSim'] = self.nSim
			self.introduction_dates['nSeq'] = self.nSeq
			self.introduction_dates['protocol'] = self.protocol

		return self.lineage_migration, self.introduction_dates

	