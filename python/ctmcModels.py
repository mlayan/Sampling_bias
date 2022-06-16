#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
MODULE CONTAINING CTMC MODELS TO SIMULATE
THE EVOLUTION OF SEQUENCES OR DISCRETE TRAITS
"""

# Description  
__author__ = 'Maylis Layan'
__creation_date__ = '2020-03-05' 
__last_update__ = '2020-03-06'

# Import libraries
import os
import sys
import math
import re
import scipy.linalg
import numpy as np
import pandas as pd
from dendropy import DnaCharacterMatrix




def sequence_sim(sequence1, timePeriod, S1, S2, eigenS):
    """
    Function to generate a new sequence according to an HKY model
        - sequence1 (str): ancestral sequence 
        - timePeriod (float): duration during which the ancestral sequence evolves
        - S1 (numpy matrix): transformation matrix of the CTMC Q matrix
        - S2 (numpy matrix): inverted S1
        - eigenS (numpy vector): eigenvalues of the CTMC Q matrix

    This function returns the new sequence as a string (same format as sequence1)
    """

    if type(sequence1) != str:
        raise TypeError('The sequence should be a string. Type = ' + type(sequence1))

    if type(timePeriod) != float and type(timePeriod) != int:
        raise TypeError('The time interval should be a float. Type = ' + type(timePeriod))

    # In the case of no elapsed time between one case and its ancestor
    if timePeriod == 0:
        return(sequence1)    

    # Compute the exponential of eigensM
    expEigens = scipy.linalg.expm(eigenS * timePeriod)
    probs = S1.dot(expEigens).dot(S2)

    transitionProbCheck = [1-1e-10] * 4
    if [sum(probs[x]) for x in range(4)] < transitionProbCheck:
        raise ValueError("Transition probabilities do not sum to 1!")

    # Proceed to the evolution of sequence1
    sequence2 = ''

    # Dictionary of correspondances
    bases = {
    'A' : 0,
    'T' : 1,
    'C' : 2,
    'G' : 3
    }

    # Iterate over the sequence
    for s in range(len(sequence1)):
        site1 = sequence1[s]
        newSite = np.random.choice(['A', 'T', 'C', 'G'], 
                                    1, 
                                    replace = False, 
                                    p = probs[bases[site1]])

        sequence2 += ''.join(newSite)

    return(sequence2)
    




def location_sim(state1, timePeriod, stateNames, 
    L1, L2, eigenL):
    """Function to compute a new state according to an asymmetric CTMC process
        - state1 (str): ancestral state 
        - timePeriod (float): duration during which the ancestral sequence evolves
        - stateNames (list of str): list of state names 
        - L1 (numpy matrix): transformation matrix of the CTMC Q matrix
        - L2 (numpy matrix): inverted L1
        - eigenL (numpy vector): eigenvalues of the CTMC Q matrix
    """

    if type(state1) != str:
        raise ValueError('The sequence should be a string. Type = ' + type(state1))

    if type(timePeriod) != float and type(timePeriod) != int:
        raise ValueError('The time interval should be a float. Type = ' + type(timePeriod))

    # In the case of no elapsed time between one case and its ancestor
    if timePeriod == 0:
        return(state1)

    # Compute the exponential of eigenL
    expEigens = scipy.linalg.expm(eigenL * timePeriod)
    probs = L1.dot(expEigens).dot(L2)

    # Check sums 
    transitionProbCheck = [1-1e-10] * len(stateNames)
    if [sum(probs[x]) for x in range(len(stateNames))] < transitionProbCheck:
        raise ValueError("Transition probabilities do not sum to 1!")

    # On-site dictionary 
    onSiteDictionary = dict(zip(stateNames, list(range(len(stateNames)))))

    # Proceed to the evolution of sequence1
    state2 = np.random.choice(stateNames,
        1, 
        replace = False, 
        p = probs[onSiteDictionary[state1]])
    state2 = ''.join(state2)

    return(state2)





def character_evolution(tree, rootSequence, rootState, discreteNames, 
    S1, S2, eigenS, L1, L2, eigenL, fileName, directory = None):
    """
    Function to simulate sequences and states along trees
        - tree (tree class from dendropy): phylogenetic tree 
        - rootSequence (str): root sequence
        - rootState (str): root state
        - discreteNames:
        - L1 (numpy matrix): transformation matrix of the CTMC Q matrix
        - L2 (numpy matrix): inverted S1
        - eigenL (numpy vector): eigenvalues of the CTMC Q matrix
        - S1 (numpy matrix): transformation matrix of the CTMC Q matrix
        - S2 (numpy matrix): inverted S1
        - eigenS (numpy vector): eigenvalues of the CTMC Q matrix
        - fileName (str): name of the file where sequences and states are written
        - directory (str): name of the directory where ouput files are written
            Default: current directory
            If it doesn't exist, it is created 
    
    States are writtent in a tab-separated file
    Sequences are writtent in a fasta file
    The function returns the annotated phylogeny

    """

    print(rootState)

    if not directory:
        directory = ''
    else:
        if not os.path.isdir(directory):
            os.mkdir(directory) 

        if not re.match('/$', directory):
            directory = directory + "/"

    # Initialize storage lists
    taxa = []
    sequences = []
    locations = []
    #age = []
    leaf = 1

    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            node.sequence = rootSequence
            node.location = rootState
        else:
            node.location = location_sim(node.parent_node.location, 
                node.edge.length, discreteNames, L1, L2, eigenL)
            node.sequence = sequence_sim(node.parent_node.sequence, 
                node.edge.length, S1, S2, eigenS)

        # Fill the character matrices
        if node.is_leaf():
            sequences.append(str(node.sequence))
            locations.append(str(node.location))
            #age.append(node.age)
            taxa.append(re.sub("'", "", str(node.taxon)))
            leaf += 1

    # Write the location file
    #print(taxa, age, type(taxa), type(age))
    #names = [t + "_" + str(a) for t, a in zip(taxa, age)]
    names = taxa
    traits = {'traits' : names, 'location' : locations}
    pd.DataFrame(traits).to_csv(directory + fileName + '.txt', 
                                header = True, index = False, 
                                sep = '\t', na_rep = 'NA')

    # Write the fasta file
    d = dict(zip(names, sequences))
    dna = DnaCharacterMatrix.from_dict(d)
    dna.write(path = directory + fileName + '.fasta', schema = 'fasta')

    return(tree)
