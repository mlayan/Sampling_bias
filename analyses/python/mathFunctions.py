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
import math
import numpy as np
from collections import Counter
from statistics import mean

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
	if extent == L:
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
	# Prior
	if sym:
		qk = (math.log(2) + nR - 1) / (nR*(nR-1)/2)
	else:
		qk = (math.log(2) + nR - 1) / (nR*(nR-1))

	# Prior odds
	priorqk = qk / (1 - qk)

	# Posterior odds
	if mean(x) == 1:
		posterior = mean(x) - (1/len(x))
		BF =  (posterior / (1-posterior)) / priorqk
	else:
		BF = (mean(x)/(1-mean(x))) / priorqk

	return(BF)










def ESS(trace, sampleInterval):
	"""Function to compute the ESS of trace
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
	if type(trace) != list:
		trace = list(trace)

	if min(trace) == max(trace):
		return(float('nan'))

	else:
		# sum of trace, excluding burn-in
		sum0 = 0.0    

		# Trace length
		traceLength = len(trace)

		# keep track of sums of trace(i)*trace(i_+ lag) for all lags, excluding burn-in
		squareLaggedSums = [0.0] * MAX_LAG
		autoCorrelation = [0.0] * MAX_LAG

		for i in range(traceLength):
			sum0 += trace[i]

			# calculate mean
			mean = sum0 / (i + 1)

			# calculate auto correlation for selected lag times
			sum1 = sum0
			sum2 = sum0
			for lagIndex in range(min(i + 1, MAX_LAG)):
				#squareLaggedSums[lagIndex] = squareLaggedSums[lagIndex] + trace[i - lagIndex] * trace[i]
				squareLaggedSums[lagIndex] += trace[i - lagIndex] * trace[i]
				""" The following line is the same approximation as in Tracer
				(valid since mean *(samples - lag), sum1, and sum2 are approximately the same)
				though a more accurate estimate would be
				autoCorrelation[lag] = m_fSquareLaggedSums.get(lag) - sum1 * sum2 """
				#autoCorrelation[lagIndex] = squareLaggedSums[lagIndex] - (sum1 + sum2) * mean + mean * mean * (i + 1 - lagIndex)
				#autoCorrelation[lagIndex] /= (i + 1 - lagIndex)
				autoCorrelation[lagIndex] = (squareLaggedSums[lagIndex] - (sum1 + sum2) * mean + mean * mean * (i + 1 - lagIndex)) / (i + 1 - lagIndex)
				sum1 -= trace[i - lagIndex]
				sum2 -= trace[lagIndex]


		maxLag = min(traceLength, MAX_LAG)
		integralOfACFunctionTimes2 = 0.0
		for lagIndex in range(maxLag):
			if lagIndex == 0:
				integralOfACFunctionTimes2 = autoCorrelation[0]
			elif lagIndex % 2 == 0:
				# fancy stopping criterion - see main comment in Tracer code of BEAST 1
				if autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex] > 0:
					integralOfACFunctionTimes2 += 2.0 * (autoCorrelation[lagIndex - 1] + autoCorrelation[lagIndex])
				else:
					# stop
					break

		# ACT
		ACT = sampleInterval * integralOfACFunctionTimes2 / autoCorrelation[0]

		# Standard Error Of The Mean
		#stdMean = math.sqrt(integralOfACFunctionTimes2 / traceLength)

		# ESS
		ESS = traceLength / (ACT / sampleInterval)

		#auto correlation time
		return(ESS)