#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
GENERAL FUNCTIONS TO HANDLE DIRECTORIES
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-05-25' 
__last_update__ = '2020-05-25'

# Import libraries
import os
import sys
import re










def checkDirectory(directory):
	"""
	Function which corrects the name of the directory
	to be compatible with functions writing data to files.
	"""

	if directory:
		if not re.search(r'/$', directory):
			directory += "/"
	else:
		directory = ""

	return(directory)
