#!/usr/bin/python 
# -*- coding: utf-8 -*-

"""
FUNCTION TO GENERATE BEAST1 and BEAST2 XML FILES 
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-01-16' 


# Import libraries
import numpy as np
import os
import pandas as pd
import re
import sys
from io import BytesIO
from datetime import datetime, date, timedelta
from pathlib import Path

# Import modules
from Bio import SeqIO
from Bio import Phylo
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from IPython.display import display
from lxml import etree
import pkgutil

# Tailored modules
from handleDirectories import checkDirectory


## To use regular expressions in xpath
regexpNS = 'http://exslt.org/regular-expressions'






# Get decimal dates
def dt_to_dec(dt):
    year_start = dt.year
    year_end = dt.year+1
    tot_days = date(year_end, 1,1) - date(year_start,1,1)
    return dt.year + (dt.day_of_year - dt.days_in_month) / tot_days.days 



def create_beast_xml(fileName, regionsFile, templateFile, logEvery, chainLength, 
	dirInput = None, dirOutput = None, outputName = None, prefix = None, 
	fixedTreeFile = None, fixedTreeDirectory= None, BSSVS = False, 
	markovJumps = False, equalDemes = False, adjustedBF = False, seed = None):
	
	##################################################
	## Load files
	# Directories to read from and write to files 
	dirInput = checkDirectory(dirInput)
	dirOutput = checkDirectory(dirOutput)

	# Load the region file 
	locations = pd.read_csv(dirInput + regionsFile, sep ="\t", header = 0)
	nrows = locations.shape[0]

	# Fasta file
	fasta_in = dirInput + fileName
	record = SeqIO.index(fasta_in, 'fasta')
	nseq = len(record)
	if nseq != nrows:sys.exit("The fasta file and the region file do not match !!")

	count = 0
	sequences = list()
	taxons = list()
	with open(fasta_in) as in_handle:
		for title, seq in SimpleFastaParser(in_handle):
			sequences.append(seq)
			taxons.append(title)
			count += 1

	# Create dataframe with taxon names, sequences, sampling locations and 
	# sampling dates
	seqDataFrame = pd.DataFrame(list(zip(taxons, sequences)), columns =['traits', 'fasta']) 
	seqDataFrame = pd.merge(seqDataFrame, locations, on = 'traits')
	seqDataFrame['dates'] = [re.search('(?<=_)[0-9]+[.]*[0-9]+$', traits).group(0) \
							for traits in seqDataFrame['traits']]

	# Sort dataframe by traits
	seqDataFrame.sort_values(by = ['traits'], inplace = True)
	seqDataFrame.reset_index(drop = True, inplace = True)

	# Filename
	if not prefix:
		prefix = ""

	if outputName:
		outputFileName = prefix + outputName
	else:
		outputFileName = prefix + fileName.replace(".fasta", "").replace("_sequences", "")

	##################################################
	# Create xmlBeast object
	newXML = xmlBeast(templateFile)
	
	# Mascot chains
	if 'mascot' in templateFile: 
		if 'mapper' in templateFile:
			newXML.mascot_mapper(seqDataFrame, outputFileName, chainLength, logEvery, equalDemes)
		
		elif 'glm' in templateFile:
			newXML.mascot_glm(seqDataFrame, outputFileName, dirInput, chainLength, logEvery)

		else: 
			newXML.mascot(seqDataFrame, outputFileName, chainLength, logEvery)

		# If estimation is performed on the tree topology
		if fixedTreeFile or fixedTreeDirectory:
			if not fixedTreeDirectory:
				fixedTreeDirectory = ""
			else:
				if not re.search(r'/$', fixedTreeDirectory):
					fixedTreeDirectory += "/"   
			tree = Phylo.read(fixedTreeDirectory + fixedTreeFile, 'newick')
			newXML.estimationOnFixedTree(outputFileName, tree)

	# Basta chains
	if 'basta' in templateFile:
		newXML.basta(
			seqDataFrame, outputFileName,
			chainLength, logEvery, BSSVS, equalDemes)

	# DTA chains
	if 'dta' in templateFile or 'beast1' in templateFile:
		newXML.dta(nseq, seqDataFrame, outputFileName,
			chainLength, logEvery, BSSVS, equalDemes, adjustedBF)

		## Add Markov Jumps
		if markovJumps:
			regions = sorted(seqDataFrame.regions.unique())
			newXML.addMarkovJumps(regions)

		## Estimation on fixed tree 
		if fixedTreeFile or fixedTreeDirectory:
			if not fixedTreeDirectory:
				fixedTreeDirectory = ""
			else:
				if not re.search(r'/$', fixedTreeDirectory):
					fixedTreeDirectory += "/"   
			
			tree = Phylo.read(fixedTreeDirectory + fixedTreeFile, 'newick')
			newXML.fixedTreeBeast1(tree)

	# Write XML
	newXML.write(outputFileName, dirOutput)

	# Write bash file
	bash_xml(outputFileName, dirOutput, logEvery, templateFile, chainLength, seed = seed)







def bash_xml(fileName, directory, logEvery, templateFile, chainLength, seed = None):
	"""
	Write a bash script to launch Beast jobs on the Pasteur cluster 
	For BEAST2 models, BEAST1 should be unloaded and BEAST2 loaded due to my default environment 
	"""

	# Write file 
	f = open(directory + fileName + ".sh", 'w')
	f.write("#!/bin/bash\n\n")
	f.write("#SBATCH -o " + fileName + ".log\n")
	f.write("#SBATCH -e " + fileName + ".err\n")

	if 'mascot' in templateFile :
		if '500' in fileName:
			f.write("#SBATCH -p mmmi -q mmmi\n")
		f.write("#SBATCH --mem=10000\n#SBATCH --mail-user=maylis.layan@pasteur.fr --mail-type=END\n\n")
		f.write("module unload beast/v1.10.4\n")
		f.write("module load beast/v2.6.6\n\n")
		f.write("srun beast -beagle -beagle_CPU -statefile " + fileName + ".states " + fileName + ".xml\n\n")

	elif 'basta' in templateFile:
		f.write("#SBATCH -p mmmi -q mmmi\n")
		f.write("#SBATCH --mem=10000\n#SBATCH --mail-user=maylis.layan@pasteur.fr --mail-type=END\n\n")
		f.write("module unload beast/v1.10.4\n")
		f.write("module load beast/v2.6.6\n\n")
		f.write("srun beast -beagle -beagle_CPU -statefile " + fileName + ".states " + fileName + ".xml\n\n")

	else:
		if seed:
			command = "srun java -Djava.library.path=/opt/gensoft/lib/beagle-lib/3.1.2/lib -jar ~/.beast/1.10.5/lib/beast.jar -save_every " + str(logEvery) + " -save_state " + fileName + ".states -seed " + seed + " " + fileName + ".xml\n\n"
			#command = "srun beast -save_every " + str(logEvery) + " -save_state " + fileName + ".states -seed " + seed + " " + fileName + ".xml\n\n"
		else:
			command = "srun java -Djava.library.path=/opt/gensoft/lib/beagle-lib/3.1.2/lib -jar ~/.beast/1.10.5/lib/beast.jar -save_every " + str(logEvery) + " -save_state " + fileName + ".states " + fileName + ".xml\n\n"
		f.write(command)

	f.write("exit 0")
	f.close()









class xmlBeast:
	"""
	Class of XML files used as entry files of BEAST
	"""

	def __init__(self, templateFile):

		parser = etree.XMLParser(remove_blank_text=True)
		data = pkgutil.get_data(__name__, "templates/" + templateFile)
		data = BytesIO(data)
		self.template = etree.parse(data, parser = parser)



	def mascot(self, seqDataFrame, filename, chainLength, logEvery):
		
		##################################################		
		# Root of the xml attribute
		root = self.template.getroot()

		data = root.find("data") # Get the data elements
		data.text = None

		##################################################
		## Loop to create string elements to add
		## as values in sequence, trait and typeTrait nodes
		typeTraitValue = ""
		traitValue = ""

		for seq in range(seqDataFrame.shape[0]):
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
		for param in root.xpath('run/state/stateNode/typeTrait'):
			param.set('value', typeTraitValue)

		for param in root.xpath('run/state/stateNode/trait'):
			param.set('value', traitValue)

		##################################################
		## Modify MCMC characteristics
		# Change chainLength and storeEvery attibutes of "run" node
		run = root.find('run')
		run.set('chainLength', str(chainLength))

		## Change file name 
		for node in root.iter('*'):
			attribs = node.attrib

			for key, value in attribs.items(): 
				if "filename" in value:
					newValue= value.replace("filename", filename)
					node.set(key, newValue)

				if key in ['storeEvery', 'logEvery'] : 
					node.set(key, str(logEvery))





	def mascot_mapper(self, seqDataFrame, filename, chainLength, logEvery, equalDemes):

		##################################################		
		# Root of the xml attribute
		root = self.template.getroot()

		data = root.find("data") # Get the data elements
		data.text = None

		##################################################
		## Loop to create string elements to add
		## as values in sequence, trait and typeTrait nodes
		typeTraitValue = ""
		traitValue = ""

		for seq in range(seqDataFrame.shape[0]):
			sequence = etree.Element("sequence", 
				id=format("seq_" + seqDataFrame['traits'][seq]), 
				spec="Sequence", 
				taxon=format(seqDataFrame['traits'][seq]), 
				totalcount="4", 
				value=format(seqDataFrame['fasta'][seq]))
			sequence.text = None
			data.append(sequence)

			traitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['dates'][seq] + ","
			typeTraitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['regions'][seq] + ","

		# Change value for trait and typeTrait nodes 
		for param in root.xpath('run/state/tree/trait'):
			param.set('value', traitValue)

		for param in root.xpath('run/distribution/distribution/distribution/dynamics/typeTrait'):
			param.set('value', typeTraitValue)

		##################################################
		## Model parametrization and MCMC priors
		# Equal deme sizes
		if equalDemes:
			for node in root.xpath('run/operator[@id="NeConstantBSSVSScaler.t:filename"]'):
				node.set("scaleAll", "true")
				node.set("scaleAllIndependently", "false")

		# Modify lambda parameter of Poisson prior
		nRegions = len(seqDataFrame.regions.unique())
		for node in root.xpath('run/distribution/distribution/prior/distr[@id="Poisson.0"]'):
			for child in node.getchildren():
				child.text = format(float(nRegions - 1))

		##################################################
		## Modify MCMC characteristics
		# Change chainLength and storeEvery attibutes of "run" node
		run = root.find('run')
		run.set('chainLength', str(chainLength))

		## Change file name and sampling frequency
		for node in root.iter('*'):
			attribs = node.attrib

			for key, value in attribs.items(): 
				if "filename" in value:
					newValue = value.replace("filename", filename)
					node.set(key, newValue)

				if key in ["logEvery", "storeEvery"]:
					node.set(key, str(logEvery))





	def mascot_glm(self, seqDataFrame, filename, dirInput, chainLength, logEvery):


		##################################################
		## Get incidence data 
		tc = [str(f) for f in Path(dirInput).glob("**/*") if re.search("transmission_chain", str(f))]
		if len(tc) > 1:
			raise ValueError("There are more than one transmission chain available in {0}".format(dirInput)) 

		# Get monthly incidence data
		tc = pd.read_csv(tc[0], sep = "\t")
		tc = tc[["infection.time", "regions"]][tc["infection.time"] > 0]
		tc["date"] = pd.to_datetime("1989-01-01") + tc['infection.time'].apply(lambda x: pd.Timedelta(x, unit='D'))
		tc["date"] = tc["date"].dt.to_period('m')
		tc = tc.groupby(by = ["date", "regions"]).count().reset_index()

		# Modify incidence frequencies to decimal dates 
		tc.sort_values(by = ['date'], ascending=False, inplace = True)
		tc["date"] = tc["date"].apply(lambda x: dt_to_dec(x))

		# Subset dataframe to sampled regions and incidence data over the 
		# sampling period
		sampled_regions = seqDataFrame.regions.unique()
		mrst = max(seqDataFrame.dates.astype('float'))
		tc = tc[(tc.date >= 1990) & (tc.date < mrst+1/12) & (tc.regions.isin(sampled_regions))]
		rateShifts = mrst - tc.date.unique()
		rateShifts = rateShifts.tolist()
		dateDict = dict(zip(tc.date.unique(), range(len(tc.date.unique()))))
		tc.date.replace(dateDict, inplace = True)
		tc = tc.pivot(index = "regions", columns = "date", values = "infection.time").replace({np.nan:0.1})

		# Retrieve incidence data
		incidenceData = []
		for col in tc.columns:     
			incidenceData += tc[col].tolist()

		##################################################
		## Root of the xml attribute
		root = self.template.getroot()

		data = root.find("data") # Get the data elements
		data.text = None

		## Loop to create string elements to add
		## as values in sequence, trait and typeTrait nodes
		typeTraitValue = ""
		traitValue = ""

		for seq in range(seqDataFrame.shape[0]):
			sequence = etree.Element("sequence", 
				id=format("seq_" + seqDataFrame['traits'][seq]), 
				spec="Sequence", 
				taxon=format(seqDataFrame['traits'][seq]), 
				totalcount="4", 
				value=format(seqDataFrame['fasta'][seq]))
			sequence.text = None
			data.append(sequence)

			traitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['dates'][seq] + ","
			typeTraitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['regions'][seq] + ","

		# Change value for trait and typeTrait nodes 
		for param in root.xpath('run/state/tree/trait'):
			param.set('value', traitValue)

		for param in root.xpath('run/distribution/distribution/distribution/dynamics/typeTrait'):
			param.set('value', typeTraitValue)

		##################################################
		## GLM covariates
		covList = root.findall(".//covariateList")

		## Migration rates covariates 
		nRegions = len(sampled_regions)
		dynamics = root.find(".//dynamics")
		dynamics.set("dimension", str(nRegions))

		for i in range(nRegions):
			for j in range(nRegions):
				if i != j:
					migName = "MigRate_" + str(i+1) + "_" + str(j+1)
					mig = np.zeros((nRegions, nRegions))
					mig[i,j] = 1
					mig = mig[~np.eye(mig.shape[0],dtype=bool)].reshape(mig.shape[0],-1)
					mig = mig.flatten(order='C').tolist()
					mig *= len(rateShifts)

					newPredictor = etree.Element(
						"covariates",
						id=migName,
						spec="beast.mascot.glmmodel.Covariate"
						)
					newPredictor.text = ' '.join(map(str, mig))
					covList[0].append(newPredictor)

		logTransform = etree.Element("transform", id="BooleanParameter", spec="parameter.BooleanParameter")
		logTransform.text = "false " * nRegions * (nRegions-1) 
		covList[0].append(logTransform)        

		## Effective population size covariates (incidence data)
		NePred = etree.Element("covariates", id="Ne", spec="beast.mascot.glmmodel.Covariate")
		NePred.text = " ".join(map(str, incidenceData))
		covList[1].append(NePred)

		logTransform = etree.Element("transform", id="BooleanParameter1", spec="parameter.BooleanParameter")
		logTransform.text = "true" 
		covList[1].append(logTransform)

		# Rate shifts
		rs = root.find(".//rateShifts")
		rs.text = " ".join(map(str, rateShifts))

		# Lower bound for Ne prior (upper bound for coalescent rate)
		dynamics.set("maxRate", "100")

		##################################################
		## Modify MCMC characteristics
		# Change chainLength and storeEvery attibutes of "run" node
		run = root.find('run')
		run.set('chainLength', str(chainLength))

		## Change file name and sampling frequency
		for node in root.iter('*'):
			attribs = node.attrib

			for key, value in attribs.items(): 
				if "filename" in value:
					newValue = value.replace("filename", filename)
					node.set(key, newValue)

				if key in ["logEvery", "storeEvery"]:
					node.set(key, str(logEvery))








	def basta(self, seqDataFrame, filename, chainLength, logEvery, BSSVS, equalDemes):

		##################################################
		## Initialize alignment element
		alignment = self.template.find("alignment") # Get the data elements
		alignment.text = None

		## Root of the template
		root = self.template.getroot()

		##################################################
		## Loop to create string elements to add
		## as values in sequence, trait and typeTrait nodes
		typeTraitValue = ""
		traitValue = ""

		for seq in range(seqDataFrame.shape[0]):
			sequence = etree.Element(
				"sequence", 
				taxon=format(seqDataFrame['traits'][seq]), 
				value=format(seqDataFrame['fasta'][seq])
				)
			sequence.text = None
			alignment.append(sequence)

			typeTraitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['regions'][seq] + ","
			traitValue += seqDataFrame['traits'][seq] + "=" + seqDataFrame['dates'][seq] + ","

		# Change value for trait and typeTrait nodes 
		for param in root.xpath('typeTraitSet'):
			param.set('value', typeTraitValue)

		for param in root.xpath('timeTraitSet'):
			param.set('value', traitValue)

		##################################################
		## Specify the dimension of the migration matrix
		nRegions = len(seqDataFrame.regions.unique())
		migDimension = {'popSizes': str(nRegions), 'rateMatrix': str(nRegions*(nRegions-1))}

		#     1- In the migration model
		for mig in root.findall('migrationModelVolz'):

			mig.set('nTypes', str(nRegions)) # Add the number of regions in nTypes attribute

			# If BSSVS, add the indicator variables
			if BSSVS:
				bssvs = etree.Element("rateMatrixFlags", spec="BooleanParameter", 
					value="true", dimension=str(nRegions*(nRegions-1)), id="rateMatrixFlags")
				mig.append(bssvs)

			for param in mig.getchildren():
				if param.tag in migDimension.keys():
					param.set('dimension', migDimension[param.tag])

		#    2- In tree initialization
		for mig in root.findall('run/init/migrationModelVolz'):
			mig.set('nTypes', str(nRegions)) # Add the number of regions in nTypes attribute

			for param in mig.getchildren():
				if param.tag in migDimension.keys():
					param.set('dimension', migDimension[param.tag])

		# Setting the population sizes to the same value
		if equalDemes:
			for param in root.xpath('run/operator[@id="PopSizeScaler"]'):
				param.set('scaleAll', "True")

		##################################################
		# To implement BSSVS 
		# Add several arguments
		if BSSVS:

			# Parameter priors
			for n in root.xpath('input[@spec="CompoundDistribution"]'):
				dist = etree.Element("distribution", spec='beast.math.distributions.Prior')
				x = etree.Element("x", spec='Sum', arg='@rateMatrixFlags')
				distr = etree.Element("distr", spec='Poisson')
				distr.set('lambda', str(nRegions-1))

				x.text = distr.text = None
				dist.append(x)
				dist.append(distr)

				n.append(dist)

			# State node
			for n in root.xpath('run/state'):
				stateNode = etree.Element("stateNode", idref="rateMatrixFlags")
				n.append(stateNode)

			# BSSVS operator
			for n in root.xpath('run/operator[@spec="MultiTypeTreeScaleVolz"]'):
				operator = etree.Element("operator", spec='BitFlipOperator', 
					id='bitFlipOperator', parameter='@rateMatrixFlags', weight="1")
				parent = n.getparent()
				parent.insert(parent.index(n)+1, operator)

		##################################################
		# Change chainLength and storeEvery attibutes of run node
		run = root.find('run')
		run.set('chainLength', str(chainLength))

		## Change the steps at which data is collected in files 
		for node in root.iter('*'):
			attribs = node.attrib

			for key, value in attribs.items(): 
				if key in ["logEvery", "storeEvery"]:
					node.set(key, str(logEvery))
				
				if "filename" in value:
					newValue = value.replace("filename", filename)
					node.set(key, newValue)
				




	def dta(self, nseq, seqDataFrame, filename, chainLength, logEvery, BSSVS, equalDemes, adjustedBF):

		##################################################
		## Modify the core of the XML file
		## Loop to add children nodes to the XML
		root = self.template.getroot()
		taxa = root.find("taxa") # Get the taxa elements
		taxa.text = None
		alignment = root.find("alignment") # Get the alignment elements
		alignment.text = None

		for seq in range(0, nseq):
			# Append taxon element
			taxon = etree.Element("taxon", id=format(seqDataFrame['traits'][seq]))
			taxa.append(taxon)

			# Specify taxon attributes, add sampling date and location
			taxonDate = etree.Element(
				"date", 
				value=format(seqDataFrame['dates'][seq]), 
				direction="forwards", units="years"
				)
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
		nRegions = len(seqDataFrame.regions.unique())
		root.find("generalSubstitutionModel/frequencies/frequencyModel/frequencies/parameter").set('dimension', str(nRegions))
		root.find("generalSubstitutionModel/rates/parameter").set('dimension', str(nRegions * (nRegions - 1)))
		root.find("generalSubstitutionModel/rateIndicator/parameter").set('dimension', str(nRegions * (nRegions - 1)))
		root.find('ancestralTreeLikelihood/frequencyModel/frequencies/parameter').set('dimension', str(nRegions))

		## Modify chain length and Poisson prior
		mcmc = root.find("mcmc")
		mcmc.set('chainLength', str(chainLength))
		poisson = mcmc.find('joint/prior/poissonPrior')
		poisson.set('mean', str(float(nRegions - 1)))
		poisson.set('offset', '0.0')

		##################################################
		## Randomize regions to compute an adjusted BF
		# Add the operator to randomize regions
		if adjustedBF:
			
			# Add tip swapping operator
			randomizeBF = etree.Element('tipStateSwapOperator', weight="2", uniformRandomization="true")
			randomizeBF.append(etree.Element('ancestralTreeLikelihood', idref='regions.treeLikelihood'))
			root.find("operators").append(randomizeBF)

			# Increase scaleFactor avec migration rates
			for node in root.xpath("operators/scaleOperator"):
				for child in node.getchildren():
					if child.attrib["idref"] == "regions.rates":
						node.set("scaleFactor", str(0.99))

		##################################################
		## Change the steps at which data is collected
		for node in root.iter("*"):
			attribs = node.attrib
			for key, value in attribs.items():
				if 'logEvery' in key:
					node.set(key, str(logEvery))

				if "filename" in value:
					newValue = value.replace("filename", filename)
					node.set(key, newValue)   




	


	def addMarkovJumps(self, regions):
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
		root = self.template.getroot()

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

		# Edit logTree element 
		# Add element to store the complete migration history in the trees file
		mJHistory = etree.Element('markovJumpsTreeLikelihood', idref='regions.treeLikelihood')
		root.find('mcmc/logTree').append(mJHistory)

		# Remove the element that stores the ancestral nodes states (redundant with the previous element)
		for toRemove in root.xpath('mcmc/logTree/trait[@name="regions.states"]'):
			toRemove.getparent().remove(toRemove)







	def fixedTreeBeast1(self, tree):
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
		root = self.template.getroot()

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





	def estimationOnFixedTree(self, outputFileName, tree):

		# Edit the run/init element
		for node in self.template.xpath('run/init[re:match(@id, "^RandomTree")]',namespaces={'re':regexpNS}):
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

		for node in self.template.xpath('run/operator'):
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



	def getroot(self):
		return self.template.getroot()

	def write(self, fileName, directory):

		self.template.write(
			directory + fileName + ".xml", 
			encoding="UTF-8", standalone="yes", pretty_print=True, xml_declaration=True
			)