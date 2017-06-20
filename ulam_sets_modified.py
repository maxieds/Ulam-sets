## Based on code written by Kevin Lui and Maxie D. Schmidt

from sys import argv
from oset import oset ## Source distribution: https://pypi.python.org/pypi/oset/
from copy import copy, deepcopy
import numpy as np

def infinity_norm(v):
    return max(abs(v[0]), abs(v[1]))


def two_norm_squared(v):
    return v[0]**2 + v[1]**2


def vector_sum(v, w): 
    if not isinstance(v, tuple): # handle the 1D case
         return v[0] + w[0]
    elif len(v) == 2 and len(w) == 2: 
         return (v[0]+w[0], v[1]+w[1])
    else: 
         raise ValueError("Only 1D and 2D vectors supported.")
##

def compute_ulam_set(n, init_vectors=[(1, 0), (0, 1)], norm=infinity_norm, return_indices = False):
    # ulam_set is the set of all ulam elements found so far
    ulam_set = oset(init_vectors)
    old_ulam_set = oset(ulam_set)

    # new_ulam is the set of all ulam elements found in the latest iteration
    new_ulam = oset(init_vectors)

    # pairwise sums is the set of all pairwise sums of all ulam elements found
    # so far that are not already in ulam_set
    pairwise_sums = oset([])
    set2index_dict = dict()
    for idx in range(0, len(init_vectors)): 
         set2index_dict[idx] = 0
    ##
    unique_counts = dict()
    for ic in init_vectors: 
         unique_counts[ic] = 1
    ##

    for nidx in range(1, n):
        # update pairwise sums by computing pairwise sums between new_ulam
        # elements and ulam_set elements and substracting ulam_set
        print nidx
        #print "ulam_set: ", ulam_set
        #old_ulam_set2 = [old_ulam_set[i] for i in range(0, len(old_ulam_set))]
        #new_ulam2 = [new_ulam[i] for i in range(0, len(new_ulam))]
        new_sums = []
        new_ulam = [(i, new_ulam[i]) for i in range(0, len(new_ulam))]
        if nidx == 0: # TODO: This only handles the case of two initial vectors
             new_sums = [vector_sum(init_vectors[0], init_vectors[1])]
        else: 
             for i in range(0, len(new_ulam)): 
                 (xidx, x) = new_ulam[i]
                 for yidx in range(0, len(ulam_set)): 
                     y = ulam_set[yidx]
                     if x == y: continue
                     new_sums += [vector_sum(x, y)]
                 ##
            ##
        ##
        new_ulam = []
        new_sums = list(oset(new_sums))
        for nsum in new_sums: 
             if nsum in unique_counts:
                  unique_counts[nsum] += 1
             else: 
                  unique_counts[nsum] = 1
        ##
        #print "pws: ", pairwise_sums
        pairwise_sums = pairwise_sums | oset(new_sums)
        #print "pws: ", pairwise_sums
        #print "ulam_set: ", ulam_set
        pairwise_sums = (pairwise_sums | ulam_set) - ulam_set
        #print "pws: ", pairwise_sums
        #print "ulam_set: ", ulam_set

        # remove elements that are not uniquely represented: 
        #print "new_sums: ", new_sums
        #print nidx, pairwise_sums
        #print "ulam_set: ", ulam_set
        #pairwise_sums_temp = deepcopy(pairwise_sums)
        #for i in range(0, len(pairwise_sums)): 
        #     v = pairwise_sums[i]
        #     if unique_counts[v] > 1: 
        #          pairwise_sums_temp.discard(v)
        #     else: 
        #          new_ulam += [v]
        #     ##
        ###
        #pairwise_sums = deepcopy(pairwise_sums_temp)


        # update new_ulam to be the set of pairwise sums of smallest norm
        smallest_norm = min([norm(pairwise_sums[xidx]) for xidx in range(0, len(pairwise_sums))])
        new_ulam = oset([pairwise_sums[i] for i in range(0, len(pairwise_sums)) \
                         if norm(pairwise_sums[i]) == smallest_norm])
                
        # handle the possibility of multiple vectors added at this nidx for time calculations:
        curidx = len(set2index_dict)
        for idx in range(0, len(new_ulam)): 
             set2index_dict[curidx + idx] = nidx
        ##
        
        # update ulam_set to include new_ulam
        #print "ulam_set: ", ulam_set
        old_ulam_set = oset(ulam_set)
        ulam_set = oset(old_ulam_set | new_ulam)
        #print "ulam_set: ", ulam_set
    if return_indices: 
         return [(set2index_dict[idx], ulam_set[idx]) for idx in range(0, len(ulam_set))]
    else: 
         return [ulam_set[idx] for idx in range(0, len(ulam_set))]

