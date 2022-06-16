#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MODULE TO ANALYZE SIMULATED TREES
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-04-21' 
__last_update__ = '2020-04-21'

# Root
indexCase = ["0_0"]

# RootDate
startDate = 1989.0

# Protocols
protocols = ["uniformS_2.5_150",
	"maxPerRegion_2.5_150",
	"maxPerRegionYear_2.5_150", 
	"uniformS_5_150",
	"maxPerRegion_5_150",
	"maxPerRegionYear_5_150",   
	"uniformS_10_150",
	"maxPerRegion_10_150",
	"maxPerRegionYear_10_150",
	"uniformS_20_150",         
	"maxPerRegion_20_150",      
	"maxPerRegionYear_20_150",
	"uniformS_50_150", 
	"maxPerRegion_50_150",
	"maxPerRegionYear_50_150", 
	"uniformS_2.5_500",
	"maxPerRegion_2.5_500",
	"maxPerRegionYear_2.5_500",
	"uniformS_5_500",
	"maxPerRegion_5_500",      
	"maxPerRegionYear_5_500",
	"uniformS_10_500",
	"maxPerRegion_10_500",      
	"maxPerRegionYear_10_500",
	"uniformS_20_500",
	"maxPerRegion_20_500",
	"maxPerRegionYear_20_500",
	"uniformS_50_500",
	"maxPerRegion_50_500",
	"maxPerRegionYear_50_500", 
	'uniform_150',
	'uniform_500',
	'biased_2.5_150',
	'biased_2.5_500',
	'biased_5_150',
	'biased_5_500',
	'biased_10_150',
	'biased_10_500',	
	'biased_20_150',
	'biased_20_500',
	'biased_50_150',
	'biased_50_500'
]

# sampleSize
sampleSize = [150, 500]

# Number of years 
nYears = 30
