#!/usr/bin/python

import copy
import word
import morph
import ends
import heapq

from sage.all import Tuples, Permutations

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
    
    self.paths_to_0 = None
    #this lists, for each vertex, the path (as a word in the base gens) to
    #the zero vertex
    
    if not self.compute_gens():
      self.connected = False
    else:
      self.connected = True
    
  def __repr__(self):
    return str(self.base_gens) + ' ' + str(self.base_gen_actions)
    
  def __str__(self):
    return 'Covering group with base gen action: ' + str(self.base_gens) + ' ' + str(self.base_gen_actions) + '\nFundamental group as base gen words: ' + str(self.gens_to_base_group)
  
  def compute_gens(self):
    """self.gens = our generators; self.gens_in_base_gens = our gens, in terms of base"""
    #do a breadth first search to get a spanning tree, and take the 
    #remaining edges
    self.fundamental_edges = {}
    undone_vertex_stack = [0]
    visited_vertices = []
    outgoing_edges = [[] for i in xrange(self.degree)]
    incoming_edges = [[] for i in xrange(self.degree)]
    self.paths_to_0 = ['' for i in xrange(self.degree)]
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
        if len(self.paths_to_0[current_vertex]) > 0 \
           and word.inverse(g) == self.paths_to_0[current_vertex][-1]:
           continue
           
        #here we get an old vertex, so this is a generator
        if target in visited_vertices:
          self.gens.append(alphabet[len(self.gens)])
          self.fundamental_edges[(current_vertex, g)] = self.gens[-1]
          self.fundamental_edges[(target, word.inverse(g))] = word.inverse(self.gens[-1])
          our_gen = self.paths_to_0[current_vertex] + g + word.inverse(self.paths_to_0[target])
          self.gens_to_base_group[self.gens[-1]] = our_gen
          self.gens_to_base_group[word.inverse(self.gens[-1])] = word.inverse(our_gen)
        else: #this edge is new
          if target in undone_vertex_stack: #don't do them twice
            continue
          self.paths_to_0[target] = self.paths_to_0[current_vertex] + g
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
  
  def rewrite_path_from_base_gens(self, w, start_vert=0):
    v = start_vert
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
  
  def path_end_vert(self, w):
    v = 0
    for g in w:
      if g.isupper():
        gi = self.base_gen_inds[g.swapcase()]
        v = self.base_gen_inverse_actions[gi][v]
      else:
        gi = self.base_gen_inds[g]
        v = self.base_gen_actions[gi][v]
    return v
  
  def rewrite_path_from_base_gens_get_end(self, w, start_vert=0):
    v = start_vert
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
    return (new_word, v)
  
  def lift(self, chain):
    if type(chain) == list:
      sublifts = [self.lift(x) for x in chain]
      ans = []
      for L in sublifts:
        if type(L) == list:
          ans += L
        else:
          ans.append(L)
      return ans
    vert_is_done = [False for i in xrange(self.degree)]
    loops = []
    while True:
      i=0
      while i<self.degree and vert_is_done[i]:
        i+=1
      if i == self.degree:
        break
      vert_is_done[i] = True
      cur_word, dest_vert = self.rewrite_path_from_base_gens_get_end(chain, i)
      while dest_vert != i:
        vert_is_done[dest_vert] = True
        sub_word, dest_vert = self.rewrite_path_from_base_gens_get_end(chain, dest_vert)
        cur_word += sub_word
      loops.append(cur_word)
    return (loops[0] if len(loops) == 1 else loops)
  
  
  
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
    
  def conjugation_action(self):
    """Compute a list of automorphisms on the covering group (in terms 
    of the covering group generators) corresponding to the coset 
    action.  It computes it by just going to every vertex and doing 
    all the generator words, and re-expressing them in the original covering gens"""
    pass


def cyclic_cover(base_gens, degree):
  """create a cyclic cover with the standard generators; 
  the first gen in base_gens is used as the gen of the cyclic group, and 
  it corresponds to the last gen in gens"""
  acts = [range(1,degree) + [0]] + [range(degree) for i in xrange(len(base_gens)-1)]
  G = FISubgroup(base_gens, acts)
  #now we need to remake stuff so it's right
  G.fundamental_edges = {}
  used_gens=0
  for i in xrange(degree):
    G.paths_to_0[i] = i*base_gens[0]
    for j in xrange(1, len(base_gens)):
      this_covering_gen = word.alphabet[used_gens]
      G.fundamental_edges[(i, base_gens[j])] = this_covering_gen
      G.fundamental_edges[(i, base_gens[j].swapcase())] = this_covering_gen.swapcase()
      G.gens_to_base_group[this_covering_gen] = (i*base_gens[0]) + base_gens[j] + (i*base_gens[0].swapcase())
      G.gens_to_base_group[this_covering_gen.swapcase()] = word.inverse(G.gens_to_base_group[this_covering_gen])
      used_gens+=1
  #do the covering gen
  this_covering_gen = word.alphabet[used_gens]
  G.fundamental_edges[(degree-1, base_gens[0])] = this_covering_gen
  G.fundamental_edges[(0, base_gens[0].swapcase())] = this_covering_gen.swapcase()
  G.gens_to_base_group[this_covering_gen] = degree*base_gens[0]
  G.gens_to_base_group[this_covering_gen.swapcase()] = degree*(base_gens[0].swapcase())
  return G
    

def all_covers(rank, deg, gens=None, verbose=1):
    P = Permutations(range(deg)).list()
    T = Tuples(P, rank)
    cover_rank = 1+deg*rank-deg
    if gens == None:
      base_gens = word.alphabet[:rank]
    else:
      base_gens = gens
    for t in T:
      G = FISubgroup(base_gens, t)
      if not G.connected:
        continue
      yield G




def cyclic_sort(L):
  """rotates the list so a minimal element is first"""
  min_i = min(range(len(L)), key=L.__getitem__)
  return L[min_i:] + L[:min_i]
  


def find_compatible_marked_edges(tripod_list, rank, ball_size):
  next_letters = word.next_letter_dict(rank)
  marked_edges = set()
  unmarked_edges = set()
  tripod_marked_edges = []
  #the edges in neither are unmarked and undetermined
  positions_done = set()
  made_choice = False
  positions_stack = [(0,'')]
  heapq.heapify(positions_stack) #probably unnecessary

  contradiction = None
  while True:
    #pop off the shortest undone word
    L, current_position = heapq.heappop(positions_stack)
    if L > ball_size:
      break
    
    print "Current position is: ", current_position
    
    #build the list of putative edges
    #it's 'relative to current_position':
    #('relative to id','relative to id, maybe with reverse at end',
    #  True=fixed (cannot move),[(tripod,pos_in_tripod),...])
    #and then each tripod gets a list of edges
    putative_edges = {}
    edges_for_tripods = [[None, None, None] for t in tripod_list]
    for ti,t in enumerate(tripod_list):
      for ei,e in enumerate(t):
        if e[0] in putative_edges:
          putative_edges[e[0]][-1].append( (ti,ei,0) )
          edges_for_tripods[ti][ei] = e[0]
        else:
          relid = current_position + e[0] #if it cancels, note this won't cancel
          relid_no_cancel = word.trim_trailing_cancel(relid)
          fixed = (relid_no_cancel in marked_edges)
          putative_edges[e[0]] = (relid_no_cancel, relid, fixed, [ (ti, ei, 0) ] )
          edges_for_tripods[ti][ei] = e[0]

    #now all initial edges are set up
    #find an edge to push out, either because it's forced 
    #to be unmarked, or because it participates in a contradiction
    contradiction = None
    while True:
      print "The list of putative edges is:"
      print putative_edges
      
      edge_to_push = None
      for e in putative_edges:
        #check if it's actually an unmarked edge
        if putative_edges[e][0] in unmarked_edges:
          print "Found the unmarked edge to push ", e
          edge_to_push = e
          break
      if edge_to_push==None:
        #we need to try to find a contradiction
        sorted_tripod_edges = [(tuple(cyclic_sort(eft)), i) for i,eft in enumerate(edges_for_tripods)]
        sorted_tripod_edges.sort()
        #now any reversed tripods should show up next to one another
        print "Trying to find a contradiction with the sorted tripod edge list:"
        print sorted_tripod_edges
        contradiction = None
        for i in xrange(len(sorted_tripod_edges)-1):
          if word.tripods_are_negative(sorted_tripod_edges[i][0], sorted_tripod_edges[i+1][0]):
            contradiction = ( sorted_tripod_edges[i][0], sorted_tripod_edges[i+1][0] )
            print "Found the contradiction ", contradiction
            break
        if contradiction != None:
          #get the edge to push out of the contradiction
          potential_pushouts = contradiction[0]
          print "Initial potential pushouts: ", potential_pushouts
          potential_pushouts = [e for e in potential_pushouts if putative_edges[e][0] not in marked_edges]
          print "Edge we are allowed to push", potential_pushouts
          lpp = len(potential_pushouts)
          if lpp == 1:
            edge_to_push = potential_pushouts[0]
          elif lpp > 0:
            #choose the one heading towards the id, if it exists
            made_choice = True
            for e in potential_pushouts:
              if len(putative_edges[e][0]) != len(putative_edges[e][1]):
                print "Chose the edge going towards id"
                edge_to_push = e
                break
            if edge_to_push == None:
              #if we're here, we can pick one, but there's no obvious choice
              #just pick the first one
              print "No obvious choice -- just pick the first"
              edge_to_push = potential_pushouts[0]
          
      if edge_to_push == None:
        #we cannot push any more; maybe there is a contradiction, but oh well
        print "Cannot push more"
        break
      else:
        #push out the edge until it splits, then re-evaluate
        e = edge_to_push
        #if we're pushing along it, it has to be unmarked
        unmarked_edges.add(putative_edges[e][0])
        #get the tripods involved
        tripods_involved = putative_edges[e][-1]
        print "Pushing out along the tripods ", tripods_involved
        pushed_edges = {}
        for ti, ei, elli in tripods_involved:
          #this is tripod_ind, end_ind, letter_ind
          #get the next letter
          new_e = e + tripod_list[ti][ei][elli+1]
          if new_e in pushed_edges:
            pushed_edges[new_e][-1].append( (ti, ei, elli+1) )
            edges_for_tripods[ti][ei] = new_e
          else:
            relid = word.edge_relative_to_id(current_position, new_e)
            relid_no_cancel = word.trim_trailing_cancel(relid)
            fixed = (relid_no_cancel in marked_edges)
            pushed_edges[new_e] = (relid_no_cancel, relid, fixed, [ (ti, ei, elli+1) ])
            edges_for_tripods[ti][ei] = new_e
        print "The new pushed edges: ",
        print pushed_edges
        for pe in pushed_edges:
          if pe in putative_edges:
            print "Pushed edge already putative"
            return
          putative_edges[pe] = pushed_edges[pe]
        del putative_edges[e]
    
    if contradiction != None:
      #there is a contradiction we cannot resolve
      print "Cannot resolve contradiction"
      return
    
    #lock in the edges
    print "Locking in putative edges"
    if made_choice:
      print "We did make a choice"
    else:
      print "We did not make a choice"
    print putative_edges
    for e in putative_edges:
      relid_no_cancel = putative_edges[e][0]
      if relid_no_cancel not in marked_edges:
        marked_edges.add(relid_no_cancel)
    #record the edges for each tripod
    for es in edges_for_tripods:
      tripod_marked_edges.append( [putative_edges[e][1] for e in es] )
    #modify the stack
    positions_done.add(current_position)
    lcp = len(current_position)
    for ell in next_letters[current_position[-1:]]:
      heapq.heappush(positions_stack, (lcp+1, current_position + ell))
    
  return marked_edges, unmarked_edges, tripod_marked_edges




def find_path_avoiding_marked_edges(desired_len, rank, marked_edges, unmarked_edges):
  """find a path which avoids the marked edges *and* edges that aren't 
  listed"""
  word_stack = ['']
  next_letters = word.next_letter_dict(rank)
  while len(word_stack) > 0:
    current_word = word_stack.pop()
    for ell in next_letters[current_word[-1:]]:
      nw = current_word + ell
      if nw in unmarked_edges:
        if len(nw) >= desired_len:
          return nw
        word_stack.append(nw)
  return None
  






