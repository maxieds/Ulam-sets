## Written by Maxie D. Schmidt on 2017.06.20 

from copy import copy, deepcopy

class OrderedSetLocal(object): 

     def __init__(self, init_vectors, ntime): 
          self.dataset = set([])
          self.elt2index_dict = dict()
          self.add_vectors(init_vectors, ntime)
     ## 
     
     def add_vectors(self, init_vectors, ntime): 
          self.dataset = self.dataset.union(set(init_vectors))
          for v in init_vectors: 
               self.elt2index_dict[v] = ntime
          ##
          return self
     ##
     
     def list(self): 
          return list(self.dataset)
     ##
     
     def get_indexed_list(self): 
          ilst = []
          for v in list(self.dataset): 
               entry = (self.elt2index_dict[v], v)
               ilst += [entry]
          ## 
          return sorted(ilst, key = lambda entry: entry[0])
     ##
     
     def union(self, rhsoset): 
          for v in rhsoset.list(): 
               vtime = rhsoset.elt2index_dict[v]
               self.add_vectors([v], ntime = vtime)
          ##
          return self
     ##

     def copy(self): 
          return deepcopy(self)
     ## 
     
## 
