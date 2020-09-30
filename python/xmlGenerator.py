#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
FUNCTION TO GENERATE BEAST1 and BEAST2 XML FILES 
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-01-16' 


# Import libraries
import os
import sys
import re
import pandas as pd
import numpy as np


# Import modules
from lxml import etree
from Bio import SeqIO
from Bio import Phylo
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from IPython.display import display

# Tailored modules
from handleDirectories import checkDirectory

## To use regular expressions in xpath
regexpNS = 'http://exslt.org/regular-expressions'










def xml_beast1(fileName, regionsFile, logEvery, chainLength, 
	dirInput = None, dirOutput = None, outputName = None, 
	prefix = None, fixedTreeFile = None, fixedTreeDirectory= None, 
	version = 1, substitutionModel = "hky", markovJumps = False, 
	mascotPriors = False):
	"""
	Function to generate based on :
		- Beast 1 template
		- Alignment file in fasta format
		- Tab-separated file with sequence location
	a Beast 1 XML file.

	Arguments are: 
		- fileName (str): name of the fasta file
		- regionsFile (str): name of the region file
		- logEvery (int) : step of the MCMC chain
		- chainLength (int): length of the MCMC chain
		- dirInput (str): path to the reading directory
		- dirOutput (str): path to the writing directory
		- outputName: name of the output XML file
		- prefix: prefix of the output XML file
		- version : version of the beast 1 xml template 
		- substitutionModel : either 'jc' or 'hky'
		- markovJumps (bool) : add elements linked to markov jumps to reconstruct the 
		number of migration events   
		- mascotPriors (bool): mascot priors on the indicator variable are used. Instead of
		having a Poisson prior with mean (n-1), with n the number of sampled locations, we use a 
		Poisson prior of parameter 1.0 and an offset of 1.0 

	"""
	##################################################
	## Template version
	templateFile = "template_beast1_"+ substitutionModel.lower() + "_v" + str(version) + ".xml"

	##################################################
	## Load files
	# Directories to read from and write to files 
	dirInput = checkDirectory(dirInput)
	dirOutput = checkDirectory(dirOutput)

	# Load the region file 
	regions = pd.read_csv(dirInput + regionsFile, sep ="\t", header = 0)
	nrows = regions.shape[0]
	ncols = regions.shape[1]
	nRegions = len(regions["regions"].unique()) # Number of regions

	# Fasta file
	fasta_in = dirInput + fileName
	record = SeqIO.index(fasta_in, 'fasta')
	nseq = len(record)
	if nseq != nrows:sys.exit("The fasta file and the region file do not match !!")

	count = 0
	fasta = list()
	names = list()
	with open(fasta_in) as in_handle:
		for title, seq in SimpleFastaParser(in_handle):
			fasta.append(seq)
			names.append(title)
			count += 1

	nChars = len(fasta[0])

	# Transform fasta and names lists into a data frame
	seqDataFrame = pd.DataFrame(list(zip(names, fasta)), columns =['traits', 'fasta']) 

	# Join the dataframe with sequences and the dataframe with regions
	seqDataFrame = pd.merge(seqDataFrame, regions, on = 'traits')

	# Retrieve dates from sequence names
	seqDataFrame['dates'] = [re.search('(?<=_)[0-9]+[.]*[0-9]+$', traits).group(0) \
								for traits in seqDataFrame['traits']]
	
	# Sort dataframe by traits 
	seqDataFrame.sort_values(by = ['traits'], inplace = True)
	seqDataFrame.reset_index(drop = True, inplace = True)

	# Parse XML template
	parser = etree.XMLParser(remove_blank_text=True)
	template = etree.parse(templateFile, parser=parser)
	root = template.getroot()

	# Parameters to modify
	if not prefix:
		prefix = ""

	if outputName:
		outputFileName = prefix + outputName
	else:
		outputFileName = re.sub(".fasta", "", fileName)
		outputFileName = prefix + re.sub("_sequences", "", outputFileName)

	##################################################
	## Modify the core of the XML file
	## Loop to add children nodes to the XML
	taxa = root.find("taxa") # Get the taxa elements
	taxa.text = None
	alignment = root.find("alignment") # Get the alignment elements
	alignment.text = None

	for seq in range(0, nseq):
		# Append taxon element
		taxon = etree.Element("taxon", id=format(seqDataFrame['traits'][seq]))
		taxa.append(taxon)

		# Specify taxon attributes, add sampling date and location
		taxonDate = etree.Element("date", 
			value=format(seqDataFrame['dates'][seq]), 
			direction="forwards", units="years")
		taxonLocation = etree.Element('attr', name='regions')
		taxonLocation.text = seqDataFrame['regions'][seq]
		taxon.append(taxonDate)
		taxon.append(taxonLocation)

		# Do the same for the alignment
		sequence = etree.Element("sequence")
		taxonSeq = etree.Element("taxon", idref=format(seqDataFrame['traits'][seq]))
		taxonSeq.tail = seqDataFrame['fasta'][seq]
		sequence.text = None
		sequence.append(taxonSeq)

		# Append the sequence
		alignment.append(sequence)

	# Edit beast/generalDataType to add sampling locations
	generalDataType = root.find('generalDataType')
	for x in sorted(seqDataFrame.regions.unique()):
		generalDataType.append(etree.Element("state", code=x))

	# Change region frequencies, rates, rates indicators and ancetsral locations 
	# dimensions
	root.find("generalSubstitutionModel/frequencies/frequencyModel/frequencies/parameter").set('dimension', str(nRegions))
	root.find("generalSubstitutionModel/rates/parameter").set('dimension', str(nRegions * (nRegions - 1)))
	root.find("generalSubstitutionModel/rateIndicator/parameter").set('dimension', str(nRegions * (nRegions - 1)))
	root.find('ancestralTreeLikelihood/frequencyModel/frequencies/parameter').set('dimension', str(nRegions))

	## Modify chain length and Poisson prior
	mcmc = root.find("mcmc")
	mcmc.set('chainLength', str(chainLength))
	poisson = mcmc.find('joint/prior/poissonPrior')
	if mascotPriors:
		poisson.set('mean', '1.0')
		poisson.set('offset', '1.0')
	else:
		poisson.set('mean', str(nRegions - 1) + '.0')
		poisson.set('offset', '0.0')

	##################################################     
	## Adapt according to the substitution site model
	# if substitutionModel.lower() == "hky":
	# 	template = HKYModel(template)
	if substitutionModel.lower() not in ["jc", "hky"]:
		raise ValueError("The substitution model passed is not correct. \
			It is {0} but it should be either jc or hky".format(substitutionModel))

	##################################################
	## Adapt according to the Markov Jumps
	if markovJumps:
		regions = sorted(seqDataFrame.regions.unique())
		template = addMarkovJumps(template, regions)

	##################################################
	# If estimation is performed on the tree topology
	## Fixed tree 
	if fixedTreeFile or fixedTreeDirectory:
		if not fixedTreeDirectory:
			fixedTreeDirectory = ""
		else:
			if not re.search(r'/$', fixedTreeDirectory):
				fixedTreeDirectory += "/"   
		
		tree = Phylo.read(fixedTreeDirectory + fixedTreeFile, 'newick')
		template = fixedTreeBeast1(template, tree)

	##################################################
	## Change the steps at which data is collected
	for nodes in root.findall("mcmc/"):
		attribs = nodes.attrib
	for key, value in attribs.items():
		if 'logEvery' in key:
			nodes.set(key, str(logEvery))

	##################################################
	## Change file name
	# Edit logger nodes to add output filename
	for nodes in root.iter('*'):
		attribs = nodes.attrib
		for key, value in attribs.items(): 
			if "filename" in value:
				newValue = value.replace("filename", outputFileName)
				nodes.set(key, newValue)   

	##################################################
	## Write files  
	template.write(dirOutput + outputFileName + ".xml", encoding="UTF-8", 
		standalone="yes", pretty_print=True, xml_declaration=True)

	# Write bash file
	bash_xml(outputFileName, dirOutput, logEvery, beastVersion = 1)










def addMarkovJumps(template, regions):
	"""
	This function modifies a regular beast1 template to add Markov Jumps 
	counts and rewards. This function is based on the lxml module.
	Arguments:
		- template (etree element): a beast1 template
		- regions: list of traits in the same order as in the "generalDataType" element

	Output:
		- edited template 
	"""
	# Get root of the template
	root = template.getroot()

	# List of strings corresponding to counts matrices 
	n = len(regions)
	regions = [x.replace(' ', '') for x in regions]
	childCount = {'regions.count':['1.0'] * n * n}
	rewards = etree.Element("rewards")
	for s in range(n):
		# Total counts
		childCount['regions.count'][s * n + s] = '0.0'

		# Region rewards
		r = ['0.0'] * n
		r[s] = '1.0'
		rewards.append(etree.Element("parameter", id='reward_' + regions[s], value=' '.join(r)))

		# Transition counts
		for d in range(n): 
			if s != d:
				m = ['0.0'] * n * n
				m[s * n + d] = '1.0'
				childCount['counts_' + regions[s] + '_' + regions[d]] = ' '.join(m)

	childCount['regions.count'] = ' '.join(childCount['regions.count'])      

	# Edit the ancestralTreeLikelihood node
	ancestral1 = root.find('ancestralTreeLikelihood')
	ancestral1.tag = 'markovJumpsTreeLikelihood'   # Change tag name
	#ancestral1.set('numberOfSimulants', "1")

	# Add counts and rewards matrices 
	for k, v in childCount.items():
		ancestral1.append(etree.Element('parameter', id=k, value=v))
	ancestral1.append(rewards)

	# Edit specific 'ancestralTreeLikelihood' tags 
	node1 = root.find('mcmc/joint/likelihood/ancestralTreeLikelihood')
	node1.tag = 'markovJumpsTreeLikelihood'
	node2 = root.find('mcmc/log/ancestralTreeLikelihood')
	node2.tag = 'markovJumpsTreeLikelihood'

	# Add element in logTree
	regionsCount = etree.Element('trait', name='regions.count', tag='regions.count')
	regionsCount.append(etree.Element('ancestralTreeLikelihood', idref='regions.treeLikelihood'))
	root.find('mcmc/logTree').append(regionsCount)

	"""
	# Add the complete history file 
	logTree = etree.Element('logTree', 
	fileName='filename.regions.history.trees.txt',
	logEvery=str(0),
	nexusFormat='true',
	sortTranslationTable='true')
	logTree.append(etree.Element('treeModel', idref='treeModel'))
	logTree.append(etree.Element('markovJumpsTreeLikelihood', idref='regions.treeLikelihood'))
	root.find('mcmc').append(logTree)
	"""

	return(template)











def HKYModel(template):
	"""
	Modification of BEAST 1 XML to include the HKY substitution model instead of the JC model
	Old function not used with the latest BEAST 1 XML templates.
	"""
	# Get the root of the xml file
	root = template.getroot()

	model = root.find("HKYModel")
	model.set('id', 'hky')
	model.find("kappa/parameter").set('id', 'kappa')
	model.find("kappa/parameter").set('lower', '0.0')
	model.find("kappa/parameter").set('value', '2.0')

	# Update site model
	site = root.find('siteModel/substitutionModel/HKYModel')
	site.set('idref', 'hky')

	# Add two children elements to operators
	scaleOp = etree.Element('scaleOperator', scaleFactor='0.75', weight='1')
	scaleOp.append(etree.Element('parameter', idref='kappa'))
	deltaEx = etree.Element('deltaExchange', delta='0.01', weight='1')
	deltaEx.append(etree.Element('parameter', idref='frequencies'))
	root.find('operators').insert(0, scaleOp)
	root.find('operators').insert(1, deltaEx)

	# Add two children to the mcmc priors
	logNormal = etree.Element('logNormalPrior', mu='1.0', offset='0.0', sigma='1.25')
	logNormal.append(etree.Element('parameter', idref='kappa'))
	dirichlet = etree.Element('dirichletPrior', alpha='1.0', sumsTo='1.0')
	dirichlet.append(etree.Element('parameter', idref='frequencies'))
	root.find('mcmc/joint/prior').insert(0, logNormal)
	root.find('mcmc/joint/prior').insert(1, dirichlet)

	# Add kappa and frequencies to log file
	for el in root.xpath('mcmc/log[@id="fileLog"]'):
		el.insert(7, etree.Element('parameter', idref='kappa'))
		el.insert(8, etree.Element('parameter', idref='frequencies'))    

	return(template)










def fixedTreeBeast1(template, tree):
	"""
	Modifies a Beast1 XML template to use a specific starting tree in Newick format
	and keep the topology fixed along the MCMC chain
	It modifies the following elements :
		- coalescentSimulator (it corresponds to the starting tree)
		- treeModel (treeModel to compute the tree likelihood)
		- operators (prevent any tree topology move)

	Arguments:
		- template (etree.Element): Beast1 XML template 
		- tree (Phylo): tree in Phylo format 

	Output: 
		- template edited

	"""
	root = template.getroot()

	# Modify the coalescent starting tree
	coalStartTree = root.find('coalescentSimulator')
	newickStartTree = etree.Element('newick', id="startingTree")
	newickStartTree.text = Phylo.Newick.BaseTree.Tree.format(tree, 'newick').replace("\n", "")
	root.replace(coalStartTree, newickStartTree)

	# Modify the treeModel element and idref treeModel
	# Modification of the tag of treeModel's child 
	# from coalescentTree to newick
	treeModel = root.find('treeModel')
	for nodes in treeModel.getchildren():
		if nodes.tag == "coalescentTree":
			nodes.tag = "newick"

	# Drop operators linked to tree rearrangements
	moves = ['wilsonBalding', 'narrowExchange', 'wideExchange', 'subtreeSlide']
	for m in moves:
		for nodes in root.iter(m):
			nodes.getparent().remove(nodes)

	return(template)










def xml_mascot(fileName, regionsFile, logEvery, chainLength, version, 
		dirInput = None, dirOutput = None, fixedTreeFile = None, 
		fixedTreeDirectory= None, 
		outputName = None, prefix = None, dtaPriors = False):
	"""
	Function to generate a MASCOT (Beast2) XML file based on :
		- Beast 2 template for MASCOT using tailored priors
		- Alignment file in fasta format
		- Tab-separated file with sequence location

	See the file MASCOT_TEMPLATES.txt for the detailed description of mascot 
	template templates, priors and operators. 

	Arguments are: 
		- fileName (str): name of the fasta file
		- regionsFile (str): name of the regions file
		- logEvery (int): step of the MCMC chain
		- chainLength (int): length of the MCMC chain
		- version (int): version of the mascot template
		- fixedTreeFile (str): name of the fixed tree file (in Newick)
		- fixesTreeDirectory (str): directory where trees are stored 
		- directory (str): directory where input files are and where the output XML file is written
		- outputName (str): name of the output XML file
		- prefix (str): prefix of the output XML file
		- templateFile (str): name of the template file 
		- dtaPriors (bool): dta priors on the indicator variable are used. Instead of
		having a Poisson prior with mean 1 and offset 1, we use a Poisson prior of parameter (n-1), with
		n the number of sampled locations, and no offset. 
	"""

	##################################################
	## Template version
	templateFile = "template_mascot_v" + str(version) + ".xml"

	##################################################
	## Load Files
	# Directory to read and write files
	dirInput = checkDirectory(dirInput)
	dirOutput = checkDirectory(dirOutput)

	# Load the region file 
	regions = pd.read_csv(dirInput + regionsFile, sep ="\t", header = 0)
	nrows = regions.shape[0]
	ncols = regions.shape[1]
	nRegions = len(regions["regions"].unique()) # Number of regions

	# Load the fasta file
	fasta_in = dirInput + fileName
	record = SeqIO.index(fasta_in, 'fasta')
	nseq = len(record)

	if nseq != nrows:sys.exit("The fasta file and the region file do not match !!")

	count = 0
	fasta = list()
	names = list()
	with open(fasta_in) as in_handle:
	    for title, seq in SimpleFastaParser(in_handle):
	        fasta.append(seq)
	        names.append(title)
	        count += 1

	nChars = len(fasta[0])

	# Transform fasta and names lists into a data frame
	seqDataFrame = pd.DataFrame(list(zip(names, fasta)), columns =['traits', 'fasta']) 

	# Join the dataframe with sequences and the dataframe with regions
	seqDataFrame = pd.merge(seqDataFrame, regions, on = 'traits')

	# Retrieve dates from sequence names
	seqDataFrame['dates'] = [re.search('(?<=_)[0-9]+[.]*[0-9]+$', traits).group(0) for traits in seqDataFrame['traits']]
	
	# Sort dataframe by traits 
	# Sort dataframe by traits 
	seqDataFrame.sort_values(by = ['traits'], inplace = True)
	seqDataFrame.reset_index(drop = True, inplace = True)

	# Parse XML template
	parser = etree.XMLParser(remove_blank_text=True)
	template = etree.parse(templateFile, parser=parser)
	root = template.getroot()

	data = root.find("data") # Get the data elements
	data.text = None
	
	# Parameters to modify
	if not prefix:
		prefix = ""

	if outputName:
		outputFileName = prefix + outputName
	else:
		outputFileName = re.sub(".fasta", "", fileName)
		outputFileName = prefix + re.sub("_sequences", "", outputFileName)

	poisson = root.find('run/distribution/distribution/prior/distr/parameter')
	if dtaPriors:
		poisson.text = str(nRegions - 1) + ".0"
	else:
		poisson.text = "1.0"

	##################################################
	## Loop to create string elements to add
	## as values in sequence, trait and typeTrait nodes
	typeTraitValue = ""
	traitValue = ""

	for seq in range(nrows):
		sequence = etree.Element("sequence", 
			id=format("seq_" + seqDataFrame['traits'][seq]), 
			spec="Sequence", 
			taxon=format(seqDataFrame['traits'][seq]), 
			totalcount="4", 
			value=format(seqDataFrame['fasta'][seq]))
		sequence.text = None
		data.append(sequence)

		typeTraitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['regions'][seq] + ","
		traitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['dates'][seq] + ","

	# Change value for trait and typeTrait nodes 
	for param in template.xpath('run/state/stateNode/typeTrait'):
		param.attrib['value'] = typeTraitValue

	for param in template.xpath('run/state/stateNode/trait'):
		param.attrib['value'] = traitValue

	##################################################
	## Modify lambda parameter of Poisson prior
	## for version 2 and 3
	if version in [2, 3]:
		poisson = template.findall('run/distribution/distribution/prior/distr/parameter')
		if len(poisson) == 1:
			poisson[0].text = format(float(nRegions - 1))
		else:
			sys.exit("The template file contains more than one Poisson prior for version {0}. \
				Verify the template file.".format(version))

	##################################################
	## Modify MCMC characteristics
	# Change the steps at which data is collected
	for node in root.findall('run/logger'):
		node.set('logEvery', str(logEvery))

	# Change chainLength and storeEvery attibutes of "run" node
	run = root.find('run')
	run.set('storeEvery', str(logEvery))
	run.set('chainLength', str(chainLength))

	# Change storeEvery attibute of "run/state" node
	for node in template.xpath('run/state'):
		node.set('storeEvery', str(logEvery))

	## Change file name 
	for nodes in root.iter('*'):
		attribs = nodes.attrib

		for key, value in attribs.items(): 
			if "filename" in value:
				newValue= value.replace("filename", outputFileName)
				nodes.set(key, newValue)

	##################################################
	# If estimation is performed on the tree topology
	if fixedTreeFile or fixedTreeDirectory:
		if not fixedTreeDirectory:
			fixedTreeDirectory = ""
		else:
			if not re.search(r'/$', fixedTreeDirectory):
				fixedTreeDirectory += "/"   
		tree = Phylo.read(fixedTreeDirectory + fixedTreeFile, 'newick')
		template = estimationOnFixedTree(template, outputFileName, tree)

	##################################################
	# Write file	
	template.write(dirOutput + prefix + outputFileName + ".xml", encoding="UTF-8", 
		standalone="yes", pretty_print=True, xml_declaration=True)

	# Write bash file
	bash_xml(outputFileName, dirOutput, logEvery, beastVersion = 2)










def xml_basta(fileName, regionsFile, logEvery, chainLength, version,  
	dirInput = None, dirOutput = None, BSSVS = True, equalDemes = False, 
	outputName = None, prefix = None):
	"""
	Function to generate a BASTA (Beast2) XML file based on :
		- Beast 2 template for BASTA
		- Alignment file in fasta format
		- Tab-separated file with sequence location

	Arguments are: 
		- fileName: name of the fasta file
		- regionsFile : name of the regions file
		- logEvery: step of the MCMC chain
		- chainLength: length of the MCMC chain
		- BSSVS: implementation of BSSVS
		- directory: directory where input files are and where the output XML file is written
		- equalDemes: if set to True for asymmetric matrices, only one parameter is estimated for the effective population size 
		- outputName: name of the output XML file
		- prefix: prefix of the output XML file
	"""
	##################################################
	## Template file 
	templateFile = "template_basta_v" + str(version) + ".xml"

	##################################################
	## Load Files
	# Directory to read and write files
	dirInput = checkDirectory(dirInput)
	dirOutput = checkDirectory(dirOutput)

	# Load the region file 
	regions = pd.read_csv(dirInput + regionsFile, sep ="\t", header = 0)
	nrows = regions.shape[0]
	nRegions = len(regions["regions"].unique()) # Number of regions

	# Load the fasta file
	fasta_in = dirInput + fileName
	record = SeqIO.index(fasta_in, 'fasta')
	nseq = len(record)

	if nseq != nrows:sys.exit("The fasta file and the region file do not match !!")

	count = 0
	fasta = list()
	names = list()
	with open(fasta_in) as in_handle:
		for title, seq in SimpleFastaParser(in_handle):
			fasta.append(seq)
			names.append(title)
			count += 1

	nChars = len(fasta[0])

	# Transform fasta and names lists into a data frame
	seqDataFrame = pd.DataFrame(list(zip(names, fasta)), columns =['traits', 'fasta']) 

	# Join the dataframe with sequences and the dataframe with regions
	seqDataFrame = pd.merge(seqDataFrame, regions, on = 'traits')

	# Retrieve dates from sequence names
	seqDataFrame['dates'] = [re.search('(?<=_)[0-9]+[.]*[0-9]+$', traits).group(0) for traits in seqDataFrame['traits']]

	# Sort dataframe by traits 
	seqDataFrame.sort_values(by = ['traits'], inplace = True)
	seqDataFrame.reset_index(drop = True, inplace = True)

	# Parse XML template
	parser = etree.XMLParser(remove_blank_text=True)
	template = etree.parse(templateFile, parser=parser)

	alignment = template.find("alignment") # Get the data elements
	alignment.text = None

	# Parameters to modify
	if not prefix:
		prefix = ""

	if outputName:
		outputFileName = prefix + outputName
	else:
		outputFileName = re.sub(".fasta", "", fileName)
		outputFileName = prefix + re.sub("_sequences", "", outputFileName)

	##################################################
	## Loop to create string elements to add
	## as values in sequence, trait and typeTrait nodes
	typeTraitValue = ""
	traitValue = ""

	for seq in range(nrows):
		sequence = etree.Element("sequence", 
			taxon=format(seqDataFrame['traits'][seq]), 
			value=format(seqDataFrame['fasta'][seq]))
		sequence.text = None
		alignment.append(sequence)

		typeTraitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['regions'][seq] + ","
		traitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['dates'][seq] + ","

	# Change value for trait and typeTrait nodes 
	for param in template.xpath('typeTraitSet'):
		param.attrib['value'] = typeTraitValue

	for param in template.xpath('timeTraitSet'):
		param.attrib['value'] = traitValue

	##################################################
	## Specify the dimension of the migration matrix
	migDimension = {'popSizes': str(nRegions), 'rateMatrix': str(nRegions*(nRegions-1))}

	#     1- In the migration model
	for mig in template.findall('migrationModelVolz'):

		mig.set('nTypes', str(nRegions)) # Add the number of regions in nTypes attribute

		# If BSSVS, add the indicator variables
		if BSSVS:
			bssvs = etree.Element("rateMatrixFlags", spec="BooleanParameter", 
				value="true", dimension=str(nRegions*(nRegions-1)), id="rateMatrixFlags")
			mig.append(bssvs)

		for param in mig.getchildren():
			if param.tag in migDimension.keys():
				param.attrib['dimension'] = migDimension[param.tag]

	#    2- In tree initialization
	for mig in template.findall('run/init/migrationModelVolz'):
		mig.set('nTypes', str(nRegions)) # Add the number of regions in nTypes attribute

		for param in mig.getchildren():
			if param.tag in migDimension.keys():
				param.attrib['dimension'] = migDimension[param.tag]

	# Setting the population sizes to the same value
	if equalDemes:
		for param in template.xpath('run/operator[@id="PopSizeScaler"]'):
			param.set('scaleAll', "True")

	##################################################
	# To implement BSSVS 
	# Add several arguments
	if BSSVS:

		# Parameter priors
		for n in template.xpath('input[@spec="CompoundDistribution"]'):
			dist = etree.Element("distribution", spec='beast.math.distributions.Prior')
			x = etree.Element("x", spec='Sum', arg='@rateMatrixFlags')
			distr = etree.Element("distr", spec='Poisson')
			distr.set('lambda', str(nRegions-1))

			x.text = distr.text = None
			dist.append(x)
			dist.append(distr)

			n.append(dist)

		# State node
		for n in template.xpath('run/state'):
			stateNode = etree.Element("stateNode", idref="rateMatrixFlags")
			n.append(stateNode)

		# BSSVS operator
		for n in template.xpath('run/operator[@spec="MultiTypeTreeScaleVolz"]'):
			operator = etree.Element("operator", spec='BitFlipOperator', 
				id='bitFlipOperator', parameter='@rateMatrixFlags', weight="1")
			parent = n.getparent()
			parent.insert(parent.index(n)+1, operator)

	##################################################
	# Change chainLength and storeEvery attibutes of run node
	run = template.find('run')
	run.set('storeEvery', str(logEvery))
	run.set('chainLength', str(chainLength))

	## Change the steps at which data is collected in files 
	for node in template.findall('run/logger'):
		node.set('logEvery', str(logEvery))
		
	## Change file name 
	# Edit logger nodes to add output filename
	for nodes in template.findall('run/logger'):
		attribs = nodes.attrib
		for key, value in attribs.items(): 
			if "filename" in value:
				newValue = value.replace("filename", outputFileName)
				nodes.set(key, newValue)
				
	# Write file
	if not prefix:
		prefix = ""

	template.write(dirOutput + prefix + outputFileName + ".xml", 
		encoding="UTF-8", standalone="yes", pretty_print=True, xml_declaration=True)

	# Write bash file
	bash_xml(outputFileName, dirOutput, logEvery, beastVersion = 2)












def bash_xml(fileName, directory, logEvery, beastVersion):
	""""""

	# Write file 
	f = open(directory + fileName + ".sh", 'w')
	f.write("#!/bin/bash\n\n")
	f.write("#SBATCH -o " + fileName + ".log\n")
	f.write("#SBATCH -e " + fileName + ".err\n")
	# f.write("#SBATCH --exclude=tars-694,tars-695\n")
	f.write("#SBATCH -p mmmi -q mmmi\n#SBATCH --mem=8000\n#SBATCH --mail-user=maylis.layan@pasteur.fr --mail-type=END\n\n")

	if beastVersion == 2:
		f.write("module unload beast/v1.10.4\n")
		f.write("module load beast/v2.6.1\n\n")
		f.write("srun beast -beagle -beagle_CPU -statefile " + fileName + ".states " + fileName + ".xml\n\n")

	else:
		f.write("srun beast -save_every " + str(logEvery) + " -save_state " + fileName + ".states  " + fileName + ".xml\n\n")

	f.write("exit 0")
	f.close()










def estimationOnFixedTree(template, outputFileName, tree):
	"""


	"""
	# Edit the run/init element
	for node in template.xpath('run/init[re:match(@id, "^RandomTree")]',namespaces={'re':regexpNS}):
		node.set('newick', Phylo.Newick.BaseTree.Tree.format(tree, 'newick').replace("\n", ""))
		node.set('IsLabelledNewick', 'true')
		node.set('adjustTipHeights', 'true')
		node.set('id', 'NewickTree.t:' + outputFileName)
		node.set('spec', 'beast.util.TreeParser')
		for child in node.getchildren():
			node.remove(child)

	# Remove tree operators 
	# Keeps the same topology
	fixedTopology1 = ["WilsonBalding", "SubtreeSlide"]
	fixedTopology2 = ["Exchange"] 
	# Keeps the same height of internal nodes
	fixedLength1 = ["Uniform"] 
	fixedLength2 = ["UpDownOperator", "ScaleOperator"]

	for node in template.xpath('run/operator'):
		if node.attrib['spec'] in fixedTopology1 + fixedLength1:
			node.getparent().remove(node)

		elif node.attrib['spec'] in fixedTopology2 and (not 'isNarrow' in node.attrib or node.attrib['isNarrow'] in ["false"]):
			node.getparent().remove(node)

		elif node.attrib['spec'] in "ScaleOperator" and ("TreeScaler" in node.attrib['id'] or "TreeRootScaler" in node.attrib['id']):
			node.getparent().remove(node)

		elif node.attrib['spec'] in "UpDownOperator":
			for child in node.getchildren():
				if "Tree" in child.attrib['idref']:
					node.getparent().remove(node)

	return(template)
