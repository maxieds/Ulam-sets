## Based on code written by Kevin Lui and Maxie D. Schmidt

from sys import argv
from copy import copy, deepcopy
import numpy as np
from OrderedSetLocal import OrderedSetLocal as oset

round2 = lambda input: round(input, 8)

def infinity_norm(v):
    return max(abs(v[0]), abs(v[1]))
    #return min(abs(v[0]), abs(v[1]))


def two_norm_squared(v):
    return v[0]**2 + v[1]**2


def vector_sum(v, w): 
    if not isinstance(v, tuple): # handle the 1D case
         return round2(v[0] + w[0])
    elif len(v) == 2 and len(w) == 2: 
         return (round2(v[0]+w[0]), round2(v[1]+w[1]))
    else: 
         raise ValueError("Only 1D and 2D vectors supported.")
##

def compute_ulam_set_fast(n, init_vectors=[(1, 0), (0, 1)], norm=two_norm_squared):
    # ulam_set is the set of all ulam elements found so far
    for (idx, ic) in enumerate(init_vectors): 
         (x, y) = ic
         init_vectors[idx] = (round2(x), round2(y))
    ##
    ulam_set = set(init_vectors)
    old_ulam_set = ulam_set.copy()

    # new_ulam is the set of all ulam elements found in the latest iteration
    new_ulam = set(init_vectors)

    # pairwise sums is the set of all pairwise sums of all ulam elements found
    # so far that are not already in ulam_set
    pairwise_sums = set([])
    unique_counts = dict()
    for ic in init_vectors: 
         unique_counts[ic] = 1
    ##

    for nidx in range(1, n + 1):
        # update pairwise sums by computing pairwise sums between new_ulam
        # elements and ulam_set elements and substracting ulam_set
        print nidx
        new_ulam = list(new_ulam)
        new_sums = [vector_sum(x, y) for x in list(old_ulam_set) for y in new_ulam
                    if x != y]
        for (xidx, x) in enumerate(new_ulam): 
             for yidx in range(xidx + 1, len(new_ulam)): 
                  y = new_ulam[yidx]
                  new_sums += [vector_sum(x, y)]
             ##
        ##
        for nsum in new_sums: 
             if nsum in unique_counts and nidx > 1:
                  unique_counts[nsum] += 1
             else: 
                  unique_counts[nsum] = 1
        ##
        pairwise_sums = pairwise_sums.union(set(new_sums))
        pairwise_sums = pairwise_sums.difference(ulam_set)

        # remove elements that are not uniquely represented: 
        pairwise_sums_temp = pairwise_sums.copy()
        for v in pairwise_sums: 
             if unique_counts[v] > 1: 
                  pairwise_sums_temp.remove(v)
        ##
        pairwise_sums = pairwise_sums_temp.copy()

        # update new_ulam to be the set of pairwise sums of smallest norm
        smallest_norm = min([norm(x) for x in list(pairwise_sums)])
        new_ulam = set([x for x in list(pairwise_sums) if norm(x) == smallest_norm])
        
        # update ulam_set to include new_ulam
        old_ulam_set = ulam_set.copy()
        ulam_set = ulam_set.union(new_ulam)

    return ulam_set
    
def compute_ulam_set(n, init_vectors=[(1, 0), (0, 1)], norm=two_norm_squared, return_indices = False):
    
    if not return_indices: # Use the fastest implementation if we don't need to keep track of times: 
         return compute_ulam_set_fast(n, init_vectors = init_vectors, norm = norm)
    ##
    
    # ulam_set is the set of all ulam elements found so far
    for (idx, ic) in enumerate(init_vectors): 
         (x, y) = ic
         init_vectors[idx] = (round2(x), round2(y))
    ##
    print "init_vectors: ", init_vectors
    ulam_set = oset(init_vectors, ntime = 0)
    old_ulam_set = ulam_set.copy()

    # new_ulam is the set of all ulam elements found in the latest iteration
    new_ulam = oset(init_vectors, ntime = 0)

    # pairwise sums is the set of all pairwise sums of all ulam elements found
    # so far that are not already in ulam_set
    pairwise_sums = set([])
    unique_counts = dict()
    for ic in init_vectors: 
         unique_counts[ic] = 1
    ##

    for nidx in range(1, n + 1):
        # update pairwise sums by computing pairwise sums between new_ulam
        # elements and ulam_set elements and substracting ulam_set
        print nidx
        new_ulam = new_ulam.list()
        new_sums = [vector_sum(x, y) for x in old_ulam_set.list() for y in new_ulam
                    if x != y]
        for (xidx, x) in enumerate(new_ulam): 
             for yidx in range(xidx + 1, len(new_ulam)): 
                  y = new_ulam[yidx]
                  new_sums += [vector_sum(x, y)]
             ##
        ##
        for nsum in new_sums: 
             if nsum in unique_counts and nidx > 1:
                  unique_counts[nsum] += 1
             else: 
                  unique_counts[nsum] = 1
        ##
        pairwise_sums = pairwise_sums.union(set(new_sums))
        pairwise_sums = pairwise_sums.difference(set(ulam_set.list()))

        # remove elements that are not uniquely represented: 
        pairwise_sums_temp = pairwise_sums.copy()
        for v in pairwise_sums: 
             if unique_counts[v] > 1: 
                  pairwise_sums_temp.remove(v)
        ##
        pairwise_sums = pairwise_sums_temp.copy()

        # update new_ulam to be the set of pairwise sums of smallest norm
        smallest_norm = min([norm(x) for x in list(pairwise_sums)])
        new_ulam = oset([x for x in list(pairwise_sums) if norm(x) == smallest_norm], ntime = nidx)
        
        # update ulam_set to include new_ulam
        old_ulam_set = ulam_set.copy()
        ulam_set = ulam_set.union(new_ulam)

    return ulam_set.get_indexed_list()


