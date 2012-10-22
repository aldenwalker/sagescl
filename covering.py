#!/usr/bin/python

import copy
import word
import morph

alphabet = list('abcdefghijklmnopqrstuvwxyz')

def perm_inverse(p):
  lp = len(p)
  inv = [0 for i in xrange(lp)]
  for i in xrange(lp):
    inv[p[i]] = i
  return inv

class FISubgroup :
  def __init__(self, base_group_gens, gen_actions):
    self.degree = len(gen_actions[0])
    self.base_gens = copy.deepcopy(base_group_gens)  #gens of base group
    self.base_gen_inds = {}
    for i in xrange(len(self.base_gens)):
      self.base_gen_inds[self.base_gens[i]] = i
    self.base_gen_actions = copy.deepcopy(gen_actions)  #permutation action of base group gens
    self.base_gen_inverse_actions = [perm_inverse(x) for x in self.base_gen_actions]
    self.rank, self.gens, self.gens_to_base_group = None, None, None
    self.fundamental_edges = None 
    #this is a dict of tuples (vert_ind, base_gen) -> our_gen which give the non-spanning
    #tree edges which give the generators for our fundamental group
    #(determining a word in our group just means recording which 
    #of these edges is crossed as we cross them (reading left to right))
    if not self.compute_gens():
      self.connected = False
    else:
      self.connected = True
    
  def __repr__(self):
    return str(self.base_gens) + ' ' + str(self.base_gen_actions)
    
  def __str__(self):
    return self.repr()
  
  def compute_gens(self):
    """self.gens = our generators; self.gens_in_base_gens = our gens, in terms of base"""
    #do a breadth first search to get a spanning tree, and take the 
    #remaining edges
    self.fundamental_edges = {}
    undone_vertex_stack = [0]
    visited_vertices = []
    outgoing_edges = [[] for i in xrange(self.degree)]
    incoming_edges = [[] for i in xrange(self.degree)]
    paths_to_0 = ['' for i in xrange(self.degree)]
    self.gens = []
    self.gens_to_base_group = {}
    while len(undone_vertex_stack) > 0:
      #load the next vertex
      current_vertex = undone_vertex_stack.pop()
      if current_vertex in visited_vertices:
        continue
      visited_vertices.append(current_vertex)
      
      #go through the potential outgoing edges
      for i in xrange(2*len(self.base_gens)):
        if i<len(self.base_gens):
          g = self.base_gens[i]
          target = self.base_gen_actions[i][current_vertex]
        else:
          real_i = i - len(self.base_gens)
          g = word.inverse(self.base_gens[real_i])
          target = self.base_gen_inverse_actions[real_i][current_vertex]
          if target == current_vertex:
            #skip this because we will already have done it
            continue
        
        #don't backtrack
        if len(paths_to_0[current_vertex]) > 0 \
           and word.inverse(g) == paths_to_0[current_vertex][-1]:
           continue
           
        #here we get an old vertex, so this is a generator
        if target in visited_vertices:
          self.gens.append(alphabet[len(self.gens)])
          self.fundamental_edges[(current_vertex, g)] = self.gens[-1]
          self.fundamental_edges[(target, word.inverse(g))] = word.inverse(self.gens[-1])
          our_gen = paths_to_0[current_vertex] + g + word.inverse(paths_to_0[target])
          self.gens_to_base_group[self.gens[-1]] = our_gen
          self.gens_to_base_group[word.inverse(self.gens[-1])] = word.inverse(our_gen)
        else: #this edge is new
          if target in undone_vertex_stack: #don't do them twice
            continue
          paths_to_0[target] = paths_to_0[current_vertex] + g
          undone_vertex_stack = [target] + undone_vertex_stack
    #ok, we are done; however, if we note that we haven't 
    #visited everything, then we know that the covering space isn't connected,
    #so indicate this by returning False
    self.rank = len(self.gens)
    if len(visited_vertices) < self.degree:
      return False
    else:
      return True
      
  def rewrite_word_from_base_gens(self, w):
    if self.gens == None:
      self.compute_gens()
    v = 0
    new_word = ''
    for i in xrange(len(w)):
      if (v, w[i]) in self.fundamental_edges:
        new_word += self.fundamental_edges[(v,w[i])]
      if w[i].isupper():
        g = self.base_gen_inds[word.inverse(w[i])]
        v = self.base_gen_inverse_actions[g][v]
      else:
        g = self.base_gen_inds[w[i]]
        v = self.base_gen_actions[g][v]
    return new_word
  
  def contains(self, w):
    v = 0
    for i in xrange(len(w)):
      if w[i].isupper():
        g = self.base_gen_inds[word.inverse(w[i])]
        v = self.base_gen_inverse_actions[g][v]
      else:
        g = self.base_gen_inds[w[i]]
        v = self.base_gen_actions[g][v]
    return v == 0
  
  def hom_from_base_hom(self, A):
    if self.gens == None:
      self.compute_gens()
    new_rules = {}
    for g in self.gens:
      new_rules[g] = A.ap(self.gens_to_base_group[g])
      new_rules[g] = self.rewrite_word_from_base_gens(new_rules[g])
    return morph.morph(new_rules)
    
  def gens_in_base_group(self):
    return [self.gens_to_base_group[x] for x in alphabet[:self.rank]]


















