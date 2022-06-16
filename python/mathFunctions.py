#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MATH FUNCTIONS TO ANALYZE BEAST OUTPUTS
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-02-27' 
__last_update__ = '2020-05-01'

# Import libraries
import os
import re
from handleDirectories import checkDirectory

import math
import numpy as np
import pandas as pd
from collections import Counter
from statistics import mean
from ess import ESS_cython2

# Variables specific to the module
# The maximum lag to compute ESS 
MAX_LAG = 2000





def hpd(x, upper_only = False, lower_only = False, conf = 95):
	"""
	The conf%-HPD is computed from a list of values x.
	By default, the upper bound and the lower bound of the HPD are returned.
	By default, the 95%-HDP is returned, it can be modified by using the conf argument.
	Boolean arguments upper_only or lower_only enable to return only the specified bound.
	"""
	# Sort x
	L = len(x)
	S = sorted(x)

	if conf > 1 :
		conf /= 100

	# Number of values of x representing conf% of the length of x  
	extent = math.floor(L * conf)

	# If the length of x is too small, the 
	# conf% hpd interval corresponds to x
	if L==0 :
		if upper_only or lower_only:
			return(float("nan"))
		else:
			return([float("nan"), float("nan")])

	elif L==1:
		if upper_only or lower_only:
			return(S[0])
		else:
			return([S[0], S[0]])

	elif extent == L:
		if upper_only:
			return(S[L-1])
		elif lower_only:
			return(S[0])
		else:
			return([S[0], S[L-1]])
	
	else:
		hpdIndex = 0
		hpdRange = None

		# Test the length of all the intervals 
		# starting from the first value of S
		for i in range(L - extent + 1):
			lower = S[i]
			upper = S[i + extent - 1]
			currRange = upper - lower

			# Update the hpd interval and 
			# the index of the left-bound value
			if not hpdRange:
				hpdRange = currRange

			if currRange < hpdRange : 
				hpdRange = currRange
				hpdIndex = i

		if upper_only : 
			return(S[hpdIndex + extent - 1])
		elif lower_only: 
			return(S[hpdIndex])
		else : 
			return([S[hpdIndex], S[hpdIndex + extent - 1]])










def KLRootPrediction(rootLocations, nR, prob = False, file = None): 
	"""
	Compute the Kullback-Leibler divergence between prior and posterior distributions for the prediction of the root
	Arguments:
		- x (list) : 
			list of root location probabilities
			list of root location for each sampled step of the MCMC chain
		- nR (int) : number of regions present in the sample, it may differ from the total number 
		of regions in simulations
		- prob (bool) : 
			if True, x corresponds to the vector of root location probabilities  (beast2 - mascot)
			if False, x corresponds to the list of root location for each sampled step of the MCMC chain (beast1)
	"""

	# Frequency list
	if prob:
		prop = [x for x in rootLocations if x > 0]
		
	else:
		prop = [x / len(rootLocations) for x in Counter(rootLocations).values()]

	# Compute log(ProbPosterior/ProbPrior)
	# ProbPrior = 1/NumberOfLocations
	# Then, sum over all locations
	KL = sum([x * math.log(x * nR) for x in prop])

	return(KL)










def bayesFactorRates(x, nR, sym = False):
	"""
	Function to compute the Bayes Factor for estimated spatial rates
	See Lemey et al., 2009
	"""

	# Prior odds
	if sym:
		qk = (nR - 1) / (nR*(nR-1)/2) #(math.log(2) + nR - 1) / (nR*(nR-1)/2)
	else:
		qk = (nR - 1) / (nR*(nR-1)) #(math.log(2) + nR - 1) / (nR*(nR-1))

	prior_odds = qk / (1 - qk)

	# BF
	if isinstance(x, list) | isinstance(x, pd.core.series.Series): # When x is provided as a list
		if mean(x) == 1:
			posterior = 1 - (1/len(x))
			BF =  (posterior / (1-posterior)) / prior_odds
		else:
			BF = (mean(x)/(1-mean(x))) / prior_odds
	else : # When x is provided as a float
		BF = (x/(1-x)) / prior_odds

	return(BF)











def bayesFactor(x, nR, sym = False):
	"""
	Function to compute the Bayes Factor for estimated spatial rates
	See Lemey et al., 2009
	"""
	# Prior
	if sym:
		qk = (nR - 1) / (nR*(nR-1)/2) # (math.log(2) + nR - 1) / (nR*(nR-1)/2)
	else:
		qk = (nR - 1) / (nR*(nR-1)) # (math.log(2) + nR - 1) / (nR*(nR-1))

	# Prior odds
	priorqk = qk / (1 - qk)

	# Posterior odds
	BF = (np.float64(mean(x)) / (1-mean(x))) / priorqk

	return(BF)








def ESS(trace, sampleInterval):
	"""
    Function to compute the ESS of trace (wrapper of ESS_V3_cython)
    Arguments:
        - trace: any object containing a list of floats or integers that 
        can be coerced into a list
        The trace SHOULD NOT contain the burn in 

        - sampleInterval: the time between two steps in the trace

    This code is based on two functions from Beast2 (see source code on Github):
        - calcStats in beast.util.LogAnalyser (https://github.com/CompEvol/beast2/blob/master/src/beast/util/LogAnalyser.java)
        - ACT in beast.core.util.ESS (https://github.com/CompEvol/beast2/blob/master/src/beast/core/util/ESS.java)
	"""

    # Change the type of trace because this function works on lists
	if type(trace) == list:
		trace = np.array(trace)
	else:
		trace = np.array(list(trace))


	if min(trace) == max(trace):
		return(float('nan'))
	else:
		return ESS_cython2(trace, sampleInterval, MAX_LAG)
