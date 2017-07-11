## Written by Maxie D. Schmidt on 6/29/2017

from sys import argv
import pprint
from ulam_sets_modified import compute_ulam_set

ROUND_PREC = 8
round2 = lambda input: round(input, ROUND_PREC)

pp = pprint.PrettyPrinter(indent = 4)
pprint = lambda input: pp.pprint(input)

X = lambda v: v[0]
Y = lambda v: v[1]

def edistsq(v, w): 
     return (v[0] - w[0])**2 + (v[1] - w[1])**2
##

def get_distance_pairs(nsteps, init_vectors): 
     (v0, v1) = map(vector, init_vectors)
     bvec_func = lambda b, bidx: ((bidx + 1) / 2.0, bidx, edistsq((nsteps + 1 - b) * v0 + b * v1, (0, 0)))
     if nsteps % 2 == 0: 
          return [bvec_func(1, 1), bvec_func(nsteps, 3)]
     else: 
          rdist_pairs = []
          for b in range(1, nsteps + 1, 2): 
               rdist_pairs += [bvec_func(b, b)]
          ##
          return rdist_pairs
     ## 
## 

if __name__ == '__main__':
     init_vector_method = argv[1]
     init_vector_coords = map(float, argv[2:-1]) + [int(argv[-1])]
     init_vectors, nsteps = [], 0
     if init_vector_method == "userdef-1d": 
          init_vectors += [init_vector_coords[0], init_vector_coords[1]]
          nsteps = init_vector_coords[2]
     elif init_vector_method == "userdef-2d" or init_vector_method == "find-minv-indices": 
          v0 = (init_vector_coords[0], init_vector_coords[1])
          v1 = (init_vector_coords[2], init_vector_coords[3])
          init_vectors += [v0, v1]
          nsteps = init_vector_coords[4]
     elif init_vector_method == "random-1d": 
          init_vectors = generate_random_vectors(vector_1d = True, numvecs = 2)
          nsteps = int(argv[2])
     elif init_vector_method == "random-2d" or init_vector_method == "none": 
          init_vectors = generate_random_vectors(vector_1d = False, numvecs = 2)
          nsteps = int(argv[2])
     else: 
          raise ValueError("Unsupported initial vector method.")
     ##
     
     (v0, v1) = init_vectors
     step_minvs = []
     for nstep in range(7, nsteps + 1): 
          dist_pairs = get_distance_pairs(nstep, init_vectors)
          dist_pairs = sorted(dist_pairs, key = lambda (bidx, bidx2, edist): edist)
          first_two_indices = [list(dist_pairs[0])[0:2], list(dist_pairs[1])[0:2]] if nstep % 2 == 0 else dist_pairs 
          step_minvs += [(nstep, first_two_indices)]
     ##
     print "Minimum Magnitude Vector Indices: "
     pprint(step_minvs)
     
     compute_ulam_set(nsteps, init_vectors, return_indices = True)

##
