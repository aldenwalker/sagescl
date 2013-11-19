"""This module implements a fatgraph class"""


from word import *
from sage.all import *
import covering
import ends

import copy




def approx_rat(x, tol=0.0000001):
  """returns a rational approximation which is closer than tol"""
  if abs(x) < tol:
    return Integer(0)
  c = continued_fraction_list(x, partial_convergents=True, nterms=10)
  for r in c[1]:
    if abs(Integer(r[0])/r[1] - x) < tol:
      return Integer(r[0])/r[1]
  return None
    

def gen_reduce(W_in):
  """reduces a list of the form [(i1,s1), (i2,s2), ...]
  by cancelling whenever we find i1==i2 and s1 = -s2"""
  W = copy.deepcopy(W_in)
  i = 0
  while i < len(W)-1:
    i1,s1 = W[i]
    i2,s2 = W[i+1]
    if i1 == i2 and s1 == -s2:
      del W[i]
      del W[i]
      if i>0:
        i -= 1
    else:
      i += 1
  return W

def gen_word_in_letters(W_in):
  return ''.join([(alphabet[i] if d>0 else alphabet[i].swapcase()) for (i,d) in W_in])


class Edge:
  def __init__(self, v0, v1, Lf, Lb):
    self.source = v0
    self.dest = v1
    self.label_forward = Lf
    self.label_backward = Lb
    self.carries_folded_edges = None
    self.carried_by_edge = None
    #the following records the list of (index, True=positive) pairs for 
    #generators of pi_1(F).  This is used while folding to find normal 
    #generators for the kernel
    self.gen_word = []
    #this records which way points towards the origin (stays None if it's not in the spanning tree)
    self.toward_basepoint = None
    
  
  def __str__(self):
    s = '(' + str(self.source) + '->' + str(self.dest) + ', ' + self.label_forward + ', ' + self.label_backward + ')'
    if self.carries_folded_edges != None:
      s += '; folded: ' + str(self.carries_folded_edges)
    if hasattr(self, 'dead'):
      s += "*dead*"
    if self.carried_by_edge != None:
      s += ' carried by ' + str(self.carried_by_edge)
    if self.gen_word != None:
      s += ' gen word: ' + str(self.gen_word)
    if self.toward_basepoint != None:
      s += ' toward basepoint: ' + str(self.toward_basepoint)
    return s
    
  def __repr__(self):
    return 'Edge(' + ','.join(map(str, [self.source, self.dest, self.label_forward, self.label_backward])) + ')'
  

class Vertex:
  """a vertex class with ordered edges; note "true" means the edge is leaving"""
  def __init__(self, edge_list, direction_list=None):
    if direction_list == None:
      self.edges = edge_list
    else:
      if len(direction_list) > 0:
        if type(direction_list[0]) == bool:
          self.edges = [(e, direction_list[i]) for i,e in enumerate(edge_list)]
        else:
          self.edges = [(e, direction_list[i]==1) for i,e in enumerate(edge_list)]
      else:
        self.edges = []
    self.carries_folded_verts = None
    self.carried_by_vert = None
    self.is_basepoint = None
    
  def __str__(self):
    s =  str(self.edges) 
    if self.carries_folded_verts != None:
      s += '; folded: ' + str(self.carries_folded_verts)
    if hasattr(self, 'dead'):
      s += "*dead*"
    if self.carried_by_vert != None:
      s += ' carried by ' + str(self.carried_by_vert)
    if self.is_basepoint:
      s += ' basepoint!'
    return s
  
  def __repr__(self):
    return 'Vertex(' + str(self.edges) + ')'
  
  def find_edge_ind(self, e, dir) :
    return self.edges.index( (e, dir) )
  
  def num_edges(self):
    return len(self.edges)
    


def full_fiber_product(F,G):
  """Returns the fiber product of F and G as a list of connected fatgraphs.
     There is no meaning to the cyclic orders on the vertices, though (just graphs)"""
  #create the list of vertices -- just all pairs
  verts = [(x,y) for y in xrange(len(G.V)) for x in xrange(len(F.V))]
  verts_indices = dict([ (verts[i], i) for i in xrange(len(verts))])
  #there is an edge for pair of edges in F and G
  edges = [(x,y) for y in xrange(len(G.E)) for y in xrange(len(F.V))]
  E = []
  V = []
  #go through the edges and install them
  for i in xrange(len(edges)):
    pass

class Fatgraph:
  def __init__(self, verts, edges):
    self.V = [x for x in verts]
    self.E = [x for x in edges]
    self.unfolded_V = None
    self.unfolded_E = None
    self.orig_gen_words = None #these are the original words which gives the generators
    self.kernel_words = None #these are the words in the original generators which 
                             #were in the kernel (before folding got rid of them)
    self.basepoint = None #if it has a basepoint vertex

  def __repr__(self):
    return 'Fatgraph(' + str(self.V) + ', ' + str(self.E) + ')'
    
  def __str__(self):
    ans = ''
    ans += "Fatgraph with " +  str(len(self.V)) + " vertices and " + str(len(self.E)) + " edges\n"
    ans += "Vertices:\n"
    for i,v in enumerate(self.V):
      ans += str(i) + ": " + str(v) + '\n'
    ans += "Edges:\n"
    for i,e in enumerate(self.E):
      ans += str(i) + ": " + str(e) + '\n'
    if self.unfolded_V != None:
      ans += 'Unfolded vertices:\n'
      for i,v in enumerate(self.unfolded_V):
        ans += str(i) + ": " + str(v) + '\n'
    if self.unfolded_E != None:
      ans += 'Unfolded edges: \n'
      for i,e in enumerate(self.unfolded_E):
        ans += str(i) + ": " + str(e) + '\n'
    if self.orig_gen_words != None:
      ans += "Original gen words: \n"
      for gw in self.orig_gen_words:
        ans += str(gw) + '\n'
    if self.kernel_words != None:
      ans += "Words in the original kernel\n"
      for kw in self.kernel_words:
        ans += str(kw) + '\n'      
    return ans
  
  def write_file_new(self, filename):
    f = open(filename, 'w')
    f.write(str(len(self.V)) + ' ' + str(len(self.E)) + '\n')
    for v in self.V:
      f.write(str(len(v.edges)))
      for (e,d) in v.edges:
        signed_ind = (e+1 if d else -(e+1))
        f.write(' ' + str(signed_ind))
      f.write('\n')
    for e in self.E:
      f.write( str(e.source) + ' ' + str(e.dest) + ' ' + e.label_forward + ' ' + e.label_backward + '\n')
    f.close()
    
  def outgoing_labels(self, ind):
    edges = self.V[ind].edges
    ans = []
    for i,direc in edges:
      if direc:  #it's outgoing
        ans.append( self.E[i].label_forward )
      else:
        ans.append( self.E[i].label_backward )
    return ans
  
  def path_has_good_partition(self, vi, paths_in):
    """determine if the parts of the path passing through vert vi
    can be broken up into a good collection"""
    paths = copy.deepcopy(paths_in)
    if type(paths[0]) == tuple:
      paths = [paths]
    
    #find all the places where the paths pass through vertex vi
    pass_thru_edges = []
    for i,path in enumerate(paths):
      for j,(ei,d) in enumerate(path):
        dvi = (self.E[e].dest if d else self.E[e].source)
        if dvi == vi:
          pass_thru_edges.append( (i,j) )
    
    #for every pair of pass-thru-edges, determine if they can be connected
    #(can they be placed on the same level?)
    connectivity_edges = []
    for i1, (pi1, i_in_p1) in enumerate(pass_thru_edges):
      for i2, (pi2, i_in_p2) in list(enumerate(pass_thru_edges))[i1:]:
        if not self.do_paths_cross(paths, (pi1, i_in_p1), (pi2, i_in_p2)):
          connectivity_edges.append( (i1,i2) )
    
    #now we need to find a set of cliques so that 
    # 1) each ok pass_thru_edge appears once
    # 2) every other pass_thru_edge appears at least once
    ok_graph = Graph(connectivity_edges)
    cliques = ok_graph.cliques_maximal()
    
    
  
  
  
  def cleanup(self, maintain_watch_edges=False):
    """removes dead edges and vertices (it's an error if they are not actually dead)"""
    #figure out where the edges and vertices will go
    dest_edges = range(len(self.E))
    dest_verts = range(len(self.V))
    num_edges = 0
    num_verts = 0
    for i in xrange(len(self.E)):
      if hasattr(self.E[i], 'dead'):
        dest_edges[i] = -1
        continue
      dest_edges[i] = num_edges
      num_edges += 1
    
    for i in xrange(len(self.V)):
      if hasattr(self.V[i], 'dead'):
        dest_verts[i] = -1
        continue
      dest_verts[i] = num_verts
      num_verts += 1
    
    for i in xrange(len(self.E)):
      if dest_edges[i] == -1:
        continue
      self.E[i].source = dest_verts[self.E[i].source]
      self.E[i].dest = dest_verts[self.E[i].dest]
      if self.E[i].carries_folded_edges != None:
        for j in self.E[i].carries_folded_edges:
          self.unfolded_E[j].carried_by_edge = dest_edges[i]
      self.E[dest_edges[i]] = self.E[i]
      
    for i in xrange(len(self.V)):
      if dest_verts[i] == -1:
        continue
      for j in xrange(len(self.V[i].edges)):
        self.V[i].edges[j] = (dest_edges[self.V[i].edges[j][0]], self.V[i].edges[j][1])
      if self.V[i].carries_folded_verts != None:
        for j in self.V[i].carries_folded_verts:
          self.unfolded_V[j].carried_by_vert = dest_verts[i]
      self.V[dest_verts[i]] = self.V[i]
    del self.E[num_edges:]
    del self.V[num_verts:]
    
    if maintain_watch_edges:
      self.watch_edges = [[dest_edges[w] for w in W] for W in self.watch_edges]

  def comb(self): 
    """alters the .towards_basepoint fields of the edges to give a spanning tree"""
    if len(self.V) == 0:
      return
    if self.basepoint == None:
      self.basepoint = 0
      self.V[0].is_basepoint = True
      for i in xrange(1, len(self.V)):
        self.V[i].is_basepoint = False
    for e in self.E:
      e.toward_basepoint = None
    done_or_stacked_vert = [False for i in xrange(len(self.V))]
    vert_stack = [self.basepoint]
    done_or_stacked_vert[self.basepoint] = True
    while len(vert_stack) > 0:
      vi = vert_stack.pop()
      v = self.V[vi]
      for e,d in v.edges:
        ovi = (self.E[e].dest if d else self.E[e].source)
        if not done_or_stacked_vert[ovi]:
          self.E[e].toward_basepoint = not d
          done_or_stacked_vert[ovi] = True
          vert_stack.append(ovi)
    return      
  
  def comb_and_find_gen_words(self):
    """combs the graph and picks edeges to be the gen words, and 
    fills in what the gen words are, and finds a basepoint"""
    self.comb()
    self.orig_gen_words = []
    for ei, e in enumerate(self.E):
      if e.toward_basepoint == None: #if it's an edge giving a generator
        #figure out the paths
        path_to_source = self.edge_path_from_basepoint_to_vert(e.source)
        path_to_dest = self.edge_path_from_basepoint_to_vert(e.dest)
        path_from_dest = [(j,not d) for (j,d) in path_to_dest[::-1]]
        full_path = path_to_source + [(ei,True)] + path_from_dest
        gen_word_from_basepoint = ''.join([(self.E[j].label_forward if d else self.E[j].label_backward) for (j,d) in path_to_source])
        gen_word_to_basepoint = ''.join([(self.E[j].label_forward if d else self.E[j].label_backward) for (j,d) in path_from_dest])
        full_gen_word = gen_word_from_basepoint + e.label_forward + gen_word_to_basepoint
        self.orig_gen_words.append( [full_path, full_gen_word] )
        #label the edge with the gen word 
        e.gen_word = [(len(self.orig_gen_words)-1, 1)]
    return 
  
  def edge_path_from_basepoint_to_vert(self, vi_in):
    ans = []
    vi = vi_in
    while not self.V[vi].is_basepoint:
      v = self.V[vi]
      for e,d in v.edges:
        if self.E[e].toward_basepoint != d:
          continue
        ans.append( (e,d) )
        vi = (self.E[e].dest if d else self.E[e].source)
        break
    return [(e, not d) for (e,d) in ans[::-1]]
  
  def all_edge_paths_from_basepoint_to_verts(self):
    """returns a list of the edge paths to the basepoint for every vertex
    (this is faster than calling the above function for every vertex"""
    edge_paths = [None for i in xrange(len(self.V))]
    edge_paths[self.basepoint] = []
    for v, vi in enumerate(self.V):
      if edge_paths[vi] != None:
        continue
      ans = []
      cur_vi = vi
      visted_vi = [vi]
      while not self.V[cur_vi].is_basepoint:
        cur_v = self.V[cur_vi]
        for e,d in cur_v.edges:
          if self.E[e].toward_basepoint != d:
            continue
          ans.append( (e,d) )
          cur_vi = (self.E[e].dest if d else self.E[e].source)
          visited_vi.append(cur_vi)
          break
      ans = [(e, not d) for (e,d) in ans[::-1]]
      visted_vi = visited_vi[::-1]
      for i in xrange(len(ans)):
        edge_paths[visted_vi[i]] = ans[:i]
    return edge_paths
  
  
  def orig_gen_word_from_basepoint_to_vert(self, vi_in):
    """returns an (unreduced) list of original gens from the origin to vertex vi, 
    by following the cominb back (this requires that there is a combing)"""
    ans = []
    vi = vi_in
    while not self.V[vi].is_basepoint:
      v = self.V[vi]
      for e,d in v.edges:
        if self.E[e].toward_basepoint == d:
          w = self.E[e].gen_word
          if not d:
            w = [(g, -d) for (g,d) in w[::-1]]
          ans.extend(w)
          vi = (self.E[e].dest if d else self.E[e].source)
          break
    return [(g, -d) for (g,d) in ans[::-1]]
  
  def unfolded_edge_pair(self, v_ind, maintain_watch_edges=False, fold_kernel_edges=True):
    """returns a pair of edge indices which have the same outgoing label,
       or None if no such pair exists"""
    if not maintain_watch_edges:
      outgoing_labels = {}
      for i in xrange(len(self.V[v_ind].edges)):
        edge_ind, edge_dir = self.V[v_ind].edges[i]
        ol = (self.E[edge_ind].label_forward if edge_dir else self.E[edge_ind].label_backward)
        ov = (self.E[edge_ind].dest if edge_dir else self.E[edge_ind].source)
        if ol in outgoing_labels:
          if fold_kernel_edges or outgoing_labels[ol][1] != ov: 
            return (outgoing_labels[ol], (i,ov))
        outgoing_labels[ol] = (i,ov)
      return None
    else:
      outgoing_labels = {}
      for i in xrange(len(self.V[v_ind].edges)):
        edge_ind, edge_dir = self.V[v_ind].edges[i]
        #if edge_ind is the last edge in its watch list, then we can't 
        #fold it
        wli = self.watch_edges_dict[edge_ind]
        if self.watch_edges[wli] == [edge_ind]:
          continue
        ol = (self.E[edge_ind].label_forward if edge_dir else self.E[edge_ind].label_backward)
        ov = (self.E[edge_ind].dest if edge_dir else self.E[edge_ind].source)
        if ol in outgoing_labels:
          if fold_kernel_edges or outgoing_labels[ol][1] != ov:
            return (outgoing_labels[ol], (i,ov))
        outgoing_labels[ol] = (i,ov)
      return None
  
  def unfolded_fatgraph_edge_pair(self, v_ind):
    """returns an ordered, consecutive pair of edge indices which have the same outgoing label"""
    lve = len(self.V[v_ind].edges)
    for i,(e1,d1) in enumerate(self.V[v_ind].edges):
      (e2,d2) = self.V[v_ind].edges[(i+1)%lve]
      ol1 = (self.E[e1].label_forward if d1 else self.E[e1].label_backward)
      ol2 = (self.E[e2].label_forward if d2 else self.E[e2].label_backward)
      if ol1 == ol2:
        #we need to check that they either go to different vertices or
        #are adjacent 
        ovi1 = (self.E[e1].dest if d1 else self.E[e1].source)
        ovi2 = (self.E[e2].dest if d2 else self.E[e2].source)
        if ovi1 != ovi2:
          return ((i,ovi1), ((i+1)%lve, ovi2))
        ov_e1_ind = self.V[ovi1].edges.index( (e1, not d1) )
        ov_e2_ind = self.V[ovi1].edges.index( (e2, not d2) )
        if ov_e1_ind == (ov_e2_ind+1)%len(self.V[ovi1].edges):
          return ((i,ovi1), ((i+1)%lve, ovi2))
    return None
  
  def is_folded(self, fatgraph_folded=False):
    for i in xrange(len(self.V)):
      if (not fatgraph_folded) and self.unfolded_edge_pair(i) != None:
        return False
      elif fatgraph_folded and self.unfolded_fatgraph_edge_pair(i) != None:
        return False
    return True      
    

  



  
  def fold(self, fatgraph_fold=False, maintain_watch_edges=False, verbose=False): 
    """returns the folded version of the fatgraph, with the folded structure; 
    if fatgraph_fold=True, only do fatgraph (surface) folds
    if maintain_watch_edges!=None, it should be a list of lists of edge indices, 
    and it'll fold, but it'll refuse to fold up the last of a watch edge list"""
    #initialize the new folded fatgraph
    new_F = Fatgraph(copy.deepcopy(self.V), copy.deepcopy(self.E))
    new_F.unfolded_V = copy.deepcopy(self.V)
    new_F.unfolded_E = copy.deepcopy(self.E)
    for i in xrange(len(new_F.V)):
      new_F.V[i].carries_folded_verts = [i]
    for i in xrange(len(new_F.E)):
      new_F.E[i].carries_folded_edges = [i]
    new_F.orig_gen_words = copy.deepcopy(self.orig_gen_words)
    
    #initialize the watch list, if necessary
    if maintain_watch_edges:
      new_F.watch_edges = copy.deepcopy(self.watch_edges)
      new_F.watch_edges_dict = dict( [ (x,i) for i in xrange(len(new_F.watch_edges)) for x in new_F.watch_edges[i]] )

    while True:
      if verbose:
        print "Current fatgraph: "
        print new_F
      
      #find an unfolded vertex
      unfolded_vert = None
      unfolded_e_p = None #which particular pair of edge indices in the vert are duplicates
      for i in xrange(len(new_F.V)):
        if hasattr(new_F.V[i], 'dead'):
          continue
        if fatgraph_fold:
          unfolded_e_p = new_F.unfolded_fatgraph_edge_pair(i)
        else:
           unfolded_e_p = new_F.unfolded_edge_pair(i, maintain_watch_edges=maintain_watch_edges, fold_kernel_edges=False)
            
        if unfolded_e_p != None:
          unfolded_vert = i
          break
      if unfolded_vert == None:
        break
      
      if verbose:
        print "Found unfolded vertex ", unfolded_vert, " with edges ", unfolded_e_p
      
      #fold the edges together
      (i1,ovi1), (i2,ovi2) = unfolded_e_p
      v = new_F.V[unfolded_vert]
      ei1, d1 = v.edges[i1]
      ei2, d2 = v.edges[i2]
      
      #just a test that should be removed
      e1 = new_F.E[ei1]
      e2 = new_F.E[ei2]
      ovi1_test = (e1.dest if v.edges[i1][1] else e1.source)
      ovi2_test = (e2.dest if v.edges[i2][1] else e2.source)
      if ovi1 != ovi1_test or ovi2 != ovi2_test:
        print "Error; other vertices aren't right"
        return 1/0

      if ovi1 == ovi2 and not fatgraph_fold:
        print "Error; other vertices shouldn't be the same"
      #we need to make sure that edge 1 points to the basepoint, if
      #one of ov1 and ov2 is the basepoint
      if self.V[ovi2].is_basepoint:
        #swap edges 1 and 2
        i1, ei1, d1, ovi1,  i2, ei2, d2, ovi2, =  i2, ei2, d2, ovi2,  i1, ei1, d1, ovi1 
       
      #get the actual edges and vertices
      ov1, ov2 = new_F.V[ovi1], new_F.V[ovi2]
      e1, e2 = new_F.E[ei1], new_F.E[ei2]            

      #add the carried edges to e1
      e1.carries_folded_edges.extend(e2.carries_folded_edges)
      e2.dead = True
      
      #if we're maintaining watch edges, then we need to 
      #remove the folded edges from the watch lists
      if maintain_watch_edges:
        wli1 = new_F.watch_edges_dict[ei1]
        wli2 = new_F.watch_edges_dict[ei2]
        if ei1 in new_F.watch_edges[wli1]:
          new_F.watch_edges[wli1].remove(ei1)        
          if verbose:
            print "Removing edge ", ei1, " from watch list ", wli1
        elif verbose:
          print "Didn't need to remove ", ei1, " from watch list ", wli1
        if ei2 in new_F.watch_edges[wli2]:
          new_F.watch_edges[wli2].remove(ei2)
          if verbose:
            print "Removing edge ", ei2, " from watch list ", wli2
        elif verbose:
          print "Didn't need to remove ", ei2, " from watch list ", wli2
          
      
      if verbose:
        print "This is edges with main indices ", v.edges[i1], ' and ', v.edges[i2]
      
      if verbose:
        print "The vertices to fold together are ", ovi1, ' and ', ovi2
      
      #remove edge 2 altogether
      #note we remove the origin first, *then* the destination, and we 
      #remember the destination index.  This ensures we know where to glue 
      #the edges from vertex 1
      del v.edges[i2]
      ov2_e_ind = ov2.edges.index( (ei2, not d2) )

      #get the *outgoing* gen word along edge 1
      x = e1.gen_word
      X = [(g, -d) for (g,d) in x[::-1]]
      if not d1:
        x,X = X,x
      #and the outgoing gen word along edge 2
      y = e2.gen_word
      Y = [(g, -d) for (g,d) in y[::-1]]
      if not d2:
        y,Y = Y,y
      
      if ovi1 != ovi2:
        #get a list of the edges that are in vertex 2, and in the 
        #correct order!
        edges_from_vert_2 = ov2.edges[ov2_e_ind+1:] \
                          + ov2.edges[:ov2_e_ind]
        
        #for all these edges, make sure they point to the right place
        #and make sure their gen words are correct
        #if x,y,w are outgoing words from edge 1, edge 2, and ov2, 
        #then w should be replaced with x^{-1}yw
        #note this will correctly handle loops
        for (et, dt) in edges_from_vert_2:
          if dt:
            new_F.E[et].source = ovi1
            new_F.E[et].gen_word = gen_reduce(X + y + new_F.E[et].gen_word)
          else:
            new_F.E[et].dest = ovi1
            new_F.E[et].gen_word = gen_reduce(new_F.E[et].gen_word + Y + x)
        #if y is a loop, then x is becoming a loop, and the label shouldn't 
        #be x; it should be Xyx
        if ovi2 == unfolded_vert:
          if d2:
            e1.gen_word = gen_reduce(X + y + x)
          else:
            e1.gen_word = gen_reduce(X + Y + x)
        
        #get the index in vertex 1 
        ov1_e_ind = ov1.edges.index( (ei1, not d1) )
        
        #place the edges from vertex 2 in order, in *front* of this index
        ov1.edges = ov1.edges[:ov1_e_ind] \
                  + edges_from_vert_2      \
                  + ov1.edges[ov1_e_ind:]
        
        #extend the list of vertices that other vert 1 is carrying
        ov1.carries_folded_verts.extend(ov2.carries_folded_verts)
        
        #kill vertex 2
        ov2.dead = True
      
      else: #if ov1 == ov2, then we don't need to get rid of vertex 2 -- just delete the incoming location
        #note this can only happen in fatgraph folding, when we don't care about the kernel
        del ov2.edges[ov2_e_ind]
    

    #the graph is folded now except for kernel items, but we need to clean it up
    new_F.cleanup(maintain_watch_edges=maintain_watch_edges)

    #now the graph is folded, except for things in the kernel
    #unless we were fatgraph folding, in which case we're just done
    new_F.kernel_words = []
    if not fatgraph_fold:
      new_F.comb()
      if verbose:
        print "Now finding things in the kernel"
        print new_F
      while True:
        unfolded_vert = None
        unfolded_e_p = None
        for i in xrange(len(new_F.V)):
          if hasattr(new_F.V[i], 'dead'):
            continue
          unfolded_e_p = new_F.unfolded_edge_pair(i, maintain_watch_edges=maintain_watch_edges, fold_kernel_edges=True)
          if unfolded_e_p != None:
            unfolded_vert = i
            break
        if unfolded_vert == None:
          break
        #get the data on the edges
        (i1,ovi1), (i2,ovi2) = unfolded_e_p
        v = new_F.V[unfolded_vert]
        ei1, d1 = v.edges[i1]
        ei2, d2 = v.edges[i2]
        if ovi1 != ovi2:
          print "Error; other vertices should be the same"
        
        #one of the edges must be an edge not in the spanning true
        #make sure that edge 2 is not in the spanning tree (so we can delete it)
        if new_F.E[ei2].toward_basepoint != None:
          #swap edges 1 and 2
          i1, ei1, d1, ovi1,  i2, ei2, d2, ovi2, =  i2, ei2, d2, ovi2,  i1, ei1, d1, ovi1 

        #get the actual edges and vertices
        ov1, ov2 = new_F.V[ovi1], new_F.V[ovi2]
        e1, e2 = new_F.E[ei1], new_F.E[ei2]            

        #add the carried edges to e1
        e1.carries_folded_edges.extend(e2.carries_folded_edges)
        e2.dead = True
      
        if verbose:
          print "Found kernel element: folding ", ei1, d1, ' and ', ei2, d2

        #if we're maintaining watch edges, then we need to 
        #remove the folded edges from the watch lists
        if maintain_watch_edges:
          wli1 = new_F.watch_edges_dict[ei1]
          wli2 = new_F.watch_edges_dict[ei2]
          if ei1 in new_F.watch_edges[wli1]:
            new_F.watch_edges[wli1].remove(ei1)        
            if verbose:
              print "Removing edge ", ei1, " from watch list ", wli1
          elif verbose:
            print "Didn't need to remove ", ei1, " from watch list ", wli1
          if ei2 in new_F.watch_edges[wli2]:
            new_F.watch_edges[wli2].remove(ei2)
            if verbose:
              print "Removing edge ", ei2, " from watch list ", wli2
          elif verbose:
            print "Didn't need to remove ", ei2, " from watch list ", wli2
            
        
        #now we just need to delete edge 2
        del v.edges[i2]
        ov2_e_ind = ov2.edges.index( (ei2, not d2) )
        del ov2.edges[ov2_e_ind]

        #to figure out the word in the kernel, we get the word z from 
        #the basepoint to unfolded_vertex.  then the kernel word 
        #is zxy^{-1}z^{-1}, where x is the word away from unfolded vert along edge1, and similarly for y
        z = new_F.orig_gen_word_from_basepoint_to_vert(unfolded_vert)
        x = e1.gen_word
        if not d1:
          x = [(g, -d) for (g,d) in x[::-1]]
        Y = e2.gen_word
        if d2:
          Y = [(g, -d) for (g,d) in Y[::-1]]
        Z = [(g, -d) for (g,d) in z[::-1]]
        new_F.kernel_words.append( gen_reduce(z+x+Y+Z) )
      #the graph is totally folded now, but we need to clean up kernel edges
      new_F.cleanup(maintain_watch_edges=maintain_watch_edges)      

    return new_F
      
  
  
  def non_injective_rectangles(self):
    """for a folded fatgraph, returns the rectangles
    which can be assembled to give any loop in the kernel"""
    
    rectangles = []
    
    #go through all the edges, and look at the original edges they carry
    #any pair of original edges they carry is a legitimate rectangle
    #a rectangle is a pair ((e1, bool), (e2,bool)), where the bools
    #record whether the edge is forward or backward
    
    for e in self.E:
      for i in xrange(len(e.carries_folded_edges)):
        e1 = e.carries_folded_edges[i]
        dir1 = (self.unfolded_E[e1].label_forward == e.label_forward)
        for j in xrange(i+1, len(e.carries_folded_edges)):
          e2 = e.carries_folded_edges[j]
          dir2 = (self.unfolded_E[e2].label_forward != e.label_forward)
          rectangles.append( ((e1, dir1), (e2, dir2)) )
          rectangles.append( ((e2, not dir2), (e1, not dir1)) )
    
    return rectangles
  
  
  def non_injective_pieces(self, do_triangles=True):
    
    rectangles = []
    
    #go through all the edges, and look at the original edges they carry
    #any pair of original edges they carry is a legitimate rectangle
    #a rectangle is a pair ((e1, bool1, nbhd1), (e2,bool, nbhd2)), where the bools
    #record whether the edge is forward or backward, and the nbhds's record
    #the edges before and after
        
    for e in self.E:
      for i in xrange(len(e.carries_folded_edges)):
        e1 = e.carries_folded_edges[i]
        dir1 = (self.unfolded_E[e1].label_forward == e.label_forward)
        e1_ivert = (self.unfolded_E[e1].source if dir1 else self.unfolded_E[e1].dest)
        e1_dvert = (self.unfolded_E[e1].dest if dir1 else self.unfolded_E[e1].source)
        preceeding_edge_ops1 = [(eo, not do) for (eo,do) in self.unfolded_V[e1_ivert].edges if (eo,do) != (e1,dir1)]
        succeeding_edge_ops1 = [(eo,do) for (eo,do) in self.unfolded_V[e1_dvert].edges if (eo,do) != (e1, not dir1)]

        for j in xrange(len(e.carries_folded_edges)):
          if j == i:
            continue
          e2 = e.carries_folded_edges[j]
          dir2 = (self.unfolded_E[e2].label_forward != e.label_forward)
          e2letter = e.label_backward
          e2_ivert = (self.unfolded_E[e2].source if dir2 else self.unfolded_E[e2].dest)
          e2_dvert = (self.unfolded_E[e2].dest if dir2 else self.unfolded_E[e2].source)
          preceeding_edge_ops2 = [(eo, not do) for (eo,do) in self.unfolded_V[e2_ivert].edges if (eo,do) != (e2,dir2)]
          succeeding_edge_ops2 = [(eo, do) for (eo,do) in self.unfolded_V[e2_dvert].edges if (eo,do) != (e2, not dir2)]
          #print (j,e2,dir2,e2letter,e2_ivert, e2_dvert)
          if e1_ivert == e2_dvert and e1_dvert == e2_ivert:
            #this is a totally dummy rectangle;
            rectangles.append( ((e1, dir1, (None, None)), (e2, dir2, (None, None))) )
            
          elif e1_ivert == e2_dvert:
            #the first edge is a dummy edge; we only need to iterate over possible other letters
            for s1 in succeeding_edge_ops1:
              for p2 in preceeding_edge_ops2:
                rectangles.append( ((e1, dir1, (None, s1)), (e2, dir2, (p2, None))) )
                
          elif e1_dvert == e2_ivert:
            #the second edge is a dummy edge; we only need to do the other letters
            for p1 in preceeding_edge_ops1:
              for s2 in succeeding_edge_ops2:
                rectangles.append( ((e1, dir1, (p1, None)), (e2, dir2, (None, s2))) )

          else:
            #both edges will be nontrivial; we need to look at all of them
            for p1 in preceeding_edge_ops1:
              for s1 in succeeding_edge_ops1:
                for p2 in preceeding_edge_ops2:
                  for s2 in succeeding_edge_ops2:
                    rectangles.append( ((e1, dir1, (p1,s1)), (e2, dir2, (p2,s2))) )
    
    #an edge is ( (v1, (i1,d1)), (v2, (i2,d2)) ), where the vi are unfolded vertices, 
    #and the wi are pairs of edges which are incoming, outgoing to the vertices v1, v2
    #the inverse of edge ((v1,(i1,d1)),(v2,(i2,d2))) is ((v2,(i2,d2)),(v1,(i1,d1)))
    #note, if v1 == v2, then we can cut the tree to reduce its size, so 
    #we may assume that v1 != v2
    #we will assume that v1 < v2 for uniqueness
    #also note, v1 and v2 have to be carried by the same folded vertex
    gluing_edges = []
    ep = [ [ ((e1,not d1), (e2,d2)) for (e1, d1) in v.edges for (e2,d2) in v.edges if (e1,d1) != (e2,d2)] for v in self.unfolded_V]
    for i,v in enumerate(self.V):
      for j in v.carries_folded_verts:
        for k in v.carries_folded_verts:
          if k <= j:
            continue
          for ep1 in ep[j]:
            for ep2 in ep[k]:
              gluing_edges.append( ((j, ep1), (k,ep2)) )
    
    #this is an argument not to use C++:
    gluing_edge_index = dict([ (gluing_edges[i], i) for i in xrange(len(gluing_edges))])
    
    #now we must build the triangles
    #this is simply all triangles (triple of distinct vertices-plus-words)
    #a triangle is a triple ((v1,ep1), (v2, ep2), (v2,ep3))
    #all vertices of the triangle are unique (otherwise we could reduce the tree)
    #and all vertices must be carried by the same folded vertex
    triangles = []
    if not do_triangles:
      return rectangles, gluing_edges, triangles
    
    for v in self.V:
      for i in v.carries_folded_verts:
        for j in v.carries_folded_verts:
          if j <= i:
            continue
          for k in v.carries_folded_verts:
            if k <= i or k==j:
              continue
            for ep1 in ep[i]:
              for ep2 in ep[j]:
                for ep3 in ep[k]:
                  triangles.append( ((i,ep1),(j,ep2),(k,ep3)) )
    
    return rectangles, gluing_edges, triangles
    
  def rectangle_boundary(self, r):
    """returns the boundary of a rectangle, as a tuple of two gluing edges"""
    ((e1, dir1, ep1), (e2, dir2, ep2)) = r
    v11 = (self.unfolded_E[e1].source if dir1 else self.unfolded_E[e1].dest)
    v12 = (self.unfolded_E[e1].dest if dir1 else self.unfolded_E[e1].source)
    v21 = (self.unfolded_E[e2].source if dir2 else self.unfolded_E[e2].dest)
    v22 = (self.unfolded_E[e2].dest if dir2 else self.unfolded_E[e2].source)

    ge1 = ( (v22, ((e2, dir2), ep2[1])), (v11, (ep1[0], (e1,dir1))) )
    ge2 = ( (v12, ((e1, dir1), ep1[1])), (v21, (ep2[0], (e2,dir2))) )

    return (ge1, ge2)
  
  def rectangle_verts(self, r):
    """returns a list (potentially with dupes) of the vertices in a rectangle"""
    ((e1, dir1, ep1), (e2, dir2, ep2)) = r
    v11 = (self.unfolded_E[e1].source if dir1 else self.unfolded_E[e1].dest)
    v12 = (self.unfolded_E[e1].dest if dir1 else self.unfolded_E[e1].source)
    v21 = (self.unfolded_E[e2].source if dir2 else self.unfolded_E[e2].dest)
    v22 = (self.unfolded_E[e2].dest if dir2 else self.unfolded_E[e2].source)
    return ([v11,v22] if v11!=v22 else [v11]) + ([v12,v21] if v12!=v21 else [v12])
  
  def kernel_element(self, required_vert=None, do_triangles=True, verbose=0):
    """return an element in the kernel"""
    #we do a linear programming problem
    R, GE, T = self.non_injective_pieces(do_triangles=do_triangles)
    
    if verbose>0:
      print "Got pieces"
      sys.stdout.flush()
    
    #there's a column for every rectangle and triangle
    num_rects = len(R)
    num_triangles = len(T)
    num_cols = num_rects + len(T)
    #there's a row for every edge, plus a row to ensure that chi = 1
    num_rows = len(GE) + 1 + (1 if required_vert!=None else 0)
    
    #p = MixedIntegerLinearProgram(solver='ppl', maximization=False)
    p = MixedIntegerLinearProgram(maximization=False)
    x = p.new_variable()
    
    #the objective is L1 norm on the rectangles
    p.set_objective( p.sum([x[i] for i in xrange(num_rects)]) )
    
    ##for each edge, add the constraint
    ##notice that every triangle and rectangle cannot contain an edge more than once
    ##so it contributes only -1, 0,or 1
    #for e in GE:
      #lf = 0
      #eb = (e[1], e[0])
      #for i,r in enumerate(R):
        #rbound = self.rectangle_boundary(r)
        #if e in rbound:
          #lf += x[i]
        #elif eb in rbound:
          #lf += (-1)*x[i]
      #for i,t in enumerate(T):
        #if e in t:
          #lf += x[num_rects + i]
        #elif eb in t:
          #lf += (-1)*v[num_rects + i]
      #p.add_constraint(lf, min=0, max=0)
    
    #build a dict which tells us what the index of an edge is
    GE_index = dict([ (GE[i], i) for i in xrange(len(GE))])
    #get a linear function for every edge
    lfs = [0 for _ in GE]
    #go through the rectangles
    for i,r in enumerate(R):
      (vep11, vep12), (vep21, vep22) = self.rectangle_boundary(r)
      if vep11[1][1] != None: # make sure it's not a dummy edge
        #determine whether to add +/-:
        if vep11[0] < vep12[0]:
          lfs[ GE_index[(vep11, vep12)] ] += x[i]
        else:
          lfs[ GE_index[(vep12, vep11)] ] += (-1)*x[i]
      if vep21[1][1] != None: # make sure it's not a dummy edge
        #determine whether to add +/-:
        if vep21[0] < vep22[0]:
          lfs[ GE_index[(vep21, vep22)] ] += x[i]
        else:
          lfs[ GE_index[(vep22, vep21)] ] += (-1)*x[i]
    
    if verbose>0:
      print "Going through triangles"
      sys.stdout.flush()
    
    #go through the triangles
    for i,t in enumerate(T):
      for j in xrange(3):
        (v1, ep1), (v2, ep2) = t[j], t[(j+1)%3]
        if v1<v2:
          lfs[ GE_index[ ((v1, ep1), (v2, ep2)) ] ] += x[num_rects + i]
        else:
          lfs[ GE_index[ ((v2, ep2), (v1, ep1)) ] ] += (-1)*x[num_rects + i]
      
    if verbose>0:
      print "Done with triangles"
      sys.stdout.flush()    
    
    #now add the linear functions
    for lf in lfs:
      p.add_constraint(lf, min=0, max=0)
  
    #also add the constraint that chi=1
    #here, triangles contribute -(1/2), normal rectangles contribute 0, 
    #and rectangles with an open end constribute 1/2
    #rectangles with two open ends contribute 1
    lf = 0
    for i,r in enumerate(R):
      coef = 0
      ge1, ge2 = self.rectangle_boundary(r)
      if ge1[0][1][1] == None:
        coef += 1
      if ge2[0][1][1] == None:
        coef += 1
      if coef > 0:
        lf += coef*x[i]
    #all of the triangles contribute -1/2
    lf += p.sum([ (-1)*x[num_rects + i] for i in xrange(num_triangles)])
    p.add_constraint(lf, min=2, max=2)
    
    #if there's a vertex that we're required to hit, then add that in as 
    #a constraint
    if required_vert != None:
      v = required_vert
      lf = 0
      print "Requiring vert ", v
      for i,r in enumerate(R):
        coef = self.rectangle_verts(r).count(v)
        if coef > 0:
          lf += coef*x[i]
      p.add_constraint(lf, min=2)
      
    if verbose > 0:
      print "Done creating LP"
      sys.stdout.flush()
      if verbose > 1:
        p.show()
     
    try:
      obj_L1 = p.solve()
    except :
      print "LP Error (probably no kernel)"
      return None
    x_solution = p.get_values(x)
    
    if verbose>0:
      print "Done solving"
      sys.stdout.flush()
    
    #now we can assemble the pieces to obtain an actual tree
    #first, we clear denominators (unecessary?), taking an integer number 
    #of copies of each piece.
    #for each edge, we'll assign a permutation of the adjacent pieces
    rat_soln = [x_solution[i] for i in xrange(num_cols)]
    if not hasattr(rat_soln[0], 'denominator'):
      rat_soln = [approx_rat(x) for x in rat_soln]
    multiplier = lcm([rat_soln[i].denominator() for i in xrange(num_cols)])
    pieces = []
    piece_indices_for_edge = {}
    for i,r in enumerate(R):
      coef = multiplier*rat_soln[i]
      if coef == 0:
        continue
      rbound = self.rectangle_boundary(r)
      piece_indices_for_edge[rbound[0]] = piece_indices_for_edge.get(rbound[0], []) + [len(pieces) + j for j in xrange(coef)]
      piece_indices_for_edge[rbound[1]] = piece_indices_for_edge.get(rbound[1], []) + [len(pieces) + j for j in xrange(coef)]
      pieces.extend(coef*[('r',r)])
    num_rect_pieces = len(pieces)
    for i,t in enumerate(T):
      coef = multiplier*rat_soln[num_rects+i]
      if coef == 0:
        continue
      tbound = ( (t[0], t[1]), (t[1],t[2]), (t[2],t[0]) )
      for tb in tbound:
        piece_indices_for_edge[tb] = piece_indices_for_edge.get(tb, []) + [len(pieces) + j for j in xrange(coef)]
      pieces.extend(coef*[('t',t)])
    #the permutation is that we match up the piece_indices_for_edge for an edge and its inverse
    
    if verbose > 0:
      print "Pieces: ", pieces
    
    #now we trace out the boundary
    #we start at some location, and follow the edges
    #the result is a list of edges with signs
    have_visted_piece = [False for i in xrange(len(pieces))]
    boundaries = []
    while True:
      #find an unvisited piece 
      try:
        start_index = have_visted_piece.index(False)
      except ValueError:
        break
      edge_list = []
      boundary = ''
      cur_ind_in_piece = 0
      cur_piece_ind = start_index
      kind,p = pieces[cur_piece_ind]
      while True:
        if verbose>0:
          print "Boundary: ",boundary
          print "Edge list: ", edge_list
          print "cur_piece_ind, cur_ind_in_piece, kind, p: ", cur_piece_ind, cur_ind_in_piece, kind, p
        #read off the boundary label, and figure out what the outgoing edge is
        have_visted_piece[cur_piece_ind] = True
        if kind=='r':
          edge_list.append( p[cur_ind_in_piece][:2] )
          if edge_list[-1][1]:
            boundary += self.unfolded_E[edge_list[-1][0]].label_forward
          else:
            boundary += self.unfolded_E[edge_list[-1][0]].label_backward
          rbound = self.rectangle_boundary(p)
          outgoing_edge = rbound[1-cur_ind_in_piece]
        else:
          outgoing_edge = (p[(cur_ind_in_piece+1)%3], p[(cur_ind_in_piece+2)%3])
        
        #get the matching edge to the outgoing edge (if it's a dummy each, turn around)
        if outgoing_edge[0][1][1] == None:
          next_kind = 'r'
          next_p = p
          outgoing_edge_inv = outgoing_edge
          next_piece_ind = cur_piece_ind
        else:
          outgoing_edge_inv = (outgoing_edge[1], outgoing_edge[0])
          piece_ind_in_edge_list = piece_indices_for_edge[outgoing_edge].index(cur_piece_ind)
          next_piece_ind = piece_indices_for_edge[outgoing_edge_inv][piece_ind_in_edge_list]
          next_kind, next_p = pieces[next_piece_ind]
          
        #get what the next piece is
        if next_kind == 'r':
          if self.rectangle_boundary(next_p)[0] == outgoing_edge_inv:
            next_ind_in_piece = 0
          else:
            next_ind_in_piece = 1
        else:
          next_ind_in_piece = next_p.index( outgoing_edge_inv[0] )
        
        if next_piece_ind == start_index and next_ind_in_piece == 0: 
          #we're done
          break
        #otherwise, loop around
        kind = next_kind
        p = next_p
        cur_ind_in_piece = next_ind_in_piece
        cur_piece_ind = next_piece_ind
          
      #we've got a complete tree now
      boundaries.append( (edge_list, boundary) )
    
    return boundaries

    
    
  def cut_along_loop(self, edge_list, verbose=0):
    """returns a fatgraph which is self cut along the embedded loop 
    edge_list, which is a list of ((i,d),...), where i is an edge index and 
    d is a direction (True = forward)."""
    #for each vertex, figure out what vertices to split it into and what edges to attach
    new_F = Fatgraph(copy.deepcopy(self.V), copy.deepcopy(self.E))
    done_verts = set()
    new_edge_translation_table = {} #this records the true index of new edge (ei,0) and (ei,1)
    for e,d in edge_list:
      v_ind = (new_F.E[e].dest if d else new_F.E[e].source)
      if v_ind in done_verts:
        continue
      done_verts.add(v_ind)
      v = new_F.V[v_ind]
      #find all the edges involving this vertex
      crossing_list = []
      edges_to_crossing_list = {}
      for j,(e2,d2) in enumerate(edge_list):
        dv = (new_F.E[e2].dest if d2 else new_F.E[e2].source)
        if dv == v_ind:
          in_v_ind = v.edges.index( (e2,not d2) )
          out_v_ind = v.edges.index( edge_list[(j+1)%len(edge_list)] )
          edges_to_crossing_list[in_v_ind] = (len(crossing_list), 'in')
          edges_to_crossing_list[out_v_ind] = (len(crossing_list), 'out')
          crossing_list.append( (j, in_v_ind, out_v_ind) )
          
      if verbose > 1:
        print "Working on cut vertex ", (v_ind, v)
        print "Found crossing_list: ", crossing_list
        print "Found edges_to_crossing_list: ", edges_to_crossing_list
      #now figure out how many vertices and what edges they get etc
      new_verts = [[] for _ in xrange(len(crossing_list)+1)]
      num_verts_seen = 0
      seen_crossing_before = set()
      current_new_vert_stack = [0] # this stores what new vertices we're working on
      for j in xrange(len(v.edges)):
        if verbose>1:
          print "Dealing with edge ", j
          print "Current new_verts: ", new_verts
        if j not in edges_to_crossing_list:
          new_verts[current_new_vert_stack[-1]].append( (v.edges[j], None) )
        else: # it IS a cutting edge
          crossing_ind, d2 = edges_to_crossing_list[j]
          #print "Found its crossing_ind, d2 = ", crossing_ind, d2
          if crossing_ind in seen_crossing_before:
            first_new_vert = current_new_vert_stack[-1]
            second_new_vert = current_new_vert_stack[-2]
            del current_new_vert_stack[-1]
          else:
            first_new_vert = current_new_vert_stack[-1]
            num_verts_seen += 1
            second_new_vert = num_verts_seen
            current_new_vert_stack.append(num_verts_seen)
            seen_crossing_before.add(crossing_ind)
          #print "Found first_new_vert and second_new_vert = ", first_new_vert, second_new_vert
          if d2 == 'in':
            new_verts[first_new_vert].append( (v.edges[j], 1) )
            new_verts[second_new_vert].append( (v.edges[j], 0) )
          else:
            new_verts[first_new_vert].append( (v.edges[j], 0) )
            new_verts[second_new_vert].append( (v.edges[j], 1) )
      
      if verbose>1:
        print "Found the new verts: ", new_verts
      #now we know what edges go into what vertices
      #we need to actually make the vertices and attach the edges
      #we don't need to add a vertex for the 0th guy
      for j in xrange(len(new_verts)):
        if j==0:
          new_vert_ind = v_ind
          new_F.V[v_ind].edges = []
        else:
          new_vert_ind = len(new_F.V)
          new_F.V.append(Vertex([]))
          
        for (e2,d2), ind2 in new_verts[j]:
          if (e2,ind2) in new_edge_translation_table:
            correct_edge_ind = new_edge_translation_table[(e2,ind2)]
          else:
            if ind2 == 0 or ind2 == None:
              correct_edge_ind = e2
              new_edge_translation_table[(e2,ind2)] = e2
            else:
              correct_edge_ind = len(new_F.E)
              new_F.E.append( Edge(0,0,'','') )
              new_edge_translation_table[(e2,ind2)] = correct_edge_ind
          if d2:
            new_F.E[correct_edge_ind].source = new_vert_ind
            new_F.E[correct_edge_ind].label_forward = new_F.E[e2].label_forward
            new_F.E[correct_edge_ind].label_backward = new_F.E[e2].label_backward
          else:
            new_F.E[correct_edge_ind].dest = new_vert_ind
          new_F.V[new_vert_ind].edges.append( (correct_edge_ind, d2) )
          
    return new_F
    
  
  def cover_with_loops_embedded(self, path, verbose=0):
    """returns a fatgraph which covers self and which contains a lift of paths 
    which is embedded, also returns the lifted paths"""
    
    #for the vertices in the path, record which time it's being hit
    hit_verts_and_times = [0 for _ in xrange(len(path))]
    last_hit_time = {}
    for i,(e,d) in enumerate(path):
      dv = (self.E[e].dest if d else self.E[e].source)
      if dv in last_hit_time:
        hit_verts_and_times[i] = (dv, last_hit_time[dv]+1)
        last_hit_time[dv] += 1
      else:
        hit_verts_and_times[i] = (dv,0)
        last_hit_time[dv] = 0
    cover_deg = max([x[-1] for x in hit_verts_and_times]) + 1
    
    #make the new vertex and edge lists
    #this produces the disconnected n-fold cover
    new_verts = []
    new_edges = []
    nE = len(self.E)
    nV = len(self.V)
    for i in xrange(cover_deg):
      this_level_verts = [Vertex( [ (e+(i*nE), d) for (e,d) in v.edges] ) for v in self.V]
      new_verts.extend( this_level_verts )
      this_level_edges = [Edge( e.source+(i*nV), e.dest+(i*nV), \
                                e.label_forward, e.label_backward ) for e in self.E]
      new_edges.extend( this_level_edges )
    
    if verbose>0:
      print "Constructed the disconnected cover"
      print new_verts
      print new_edges

    #now follow the path around, transposing levels as necessary
    current_level = hit_verts_and_times[-1][-1]  #the current level needs to start on the level of the final vertex
    current_cover_edge = current_level*nE + path[0][0]
    cover_edge_path = []
    for i,(e,d) in enumerate(path):
      dv,target_level = hit_verts_and_times[i]
      cover_edge_path.append((current_cover_edge, d))
      if verbose>0:
        print "Current fatgraph: "
        print Fatgraph(new_verts, new_edges)
        print "Index ",i," in the path, edge: ", (e,d), ", with vertex target: ", (dv,target_level)
        print "Current cover edge: ", current_cover_edge

      #otherwise, we need to swap the edges around
      #this is the index in the vertex edge list 
      v_e_index = self.V[dv].edges.index( (e, not d) )
      #this is where we are currently headed
      cover_dv_ind = (new_edges[current_cover_edge].dest if d else new_edges[current_cover_edge].source)
      #this is where we *want* to be going
      cover_odv_ind = target_level*nV + dv
      #this is the index of the other covering edge
      cover_oe_ind = new_verts[cover_odv_ind].edges[ v_e_index ][0]
      
      if verbose>0:
        print "index in vert edge list: ", v_e_index
        print "current dest vertex: ", cover_dv_ind
        print "desired dest vertex: ", cover_odv_ind
        print "other edge involded: ", cover_oe_ind
      
      if cover_dv_ind != cover_odv_ind:  #maybe we're lucky and don't have to do anything, so we can skip this
        if d:
          new_edges[current_cover_edge].dest = cover_odv_ind
          new_edges[cover_oe_ind].dest = cover_dv_ind
        else:
          new_edges[current_cover_edge].source = cover_odv_ind
          new_edges[cover_oe_ind].source = cover_dv_ind
        new_verts[cover_dv_ind].edges[ v_e_index ] = (cover_oe_ind, not d)
        new_verts[cover_odv_ind].edges[ v_e_index ] = (current_cover_edge, not d)
      
      next_e,next_d = path[ (i+1)%len(path) ]
      next_outgoing_v_e_index = self.V[dv].edges.index( (next_e, next_d) )
      current_cover_edge = new_verts[cover_odv_ind].edges[ next_outgoing_v_e_index ][0]

    
    return (Fatgraph(new_verts, new_edges), cover_edge_path)
    
    
  def small_cover_with_loops_embedded(self, paths, verbose=0):
    """construct a cover of a smallish degree in which the paths are embedded; 
    it returns the covering fatgraph, plus the covering paths in the cover"""
    
    #get the collection of edge pairs around every vertex
    #an edge pair is (i,j), for the ith path, jth and j+1st position
    edge_pairs = [[] for _ in xrange(len(self.V))]
    for i,p in enumerate(paths):
      for j,(e,d) in enumerate(p):
        vi = (self.E[e].dest if d else self.E[e].source)
        edge_pairs[vi].append( (i,j) )
    
    #greedily build collections of edge pairs to appear on 
    #the various levels at each vertex
    #the groups_of_edge_pairs will be collections of sets of indices 
    #into the edge_pairs[i] list
    groups_of_edge_pairs = [[] for _ in xrange(len(self.V))]
    for i,v in enumerate(self.V):
      done_eps = [False for _ in xrange(len(edge_pairs[i]))]
      while False in done_eps:
        current_ep_set = set()
        #scan through and add all the edge pairs we can
        #this only requires one scan through
        for epi1, (j,k) in enumerate(edge_pairs[i]):
          if all([are_compatible_eps( (j,k), edge_pairs[i][epi2] ) for epi2 in current_ep_set]):
            current_ep_set.add( (j,k) )
            done_eps[epi1] = True
        #add the current set
        groups_of_edge_pairs[i].append( current_ep_set )
    
    #now we have collections of edge pairs
    cover_degree = max(map(len, groups_of_edge_pairs))
    
    #for simplicity, we'll just record the level of all the edge pairs
    edge_pair_levels = {}
    for i,v in enumerate(self.V):
      for j,gep in enumerate(groups_of_edge_pairs[i]):
        for ep in gep:
          if ep in edge_pair_levels:
            print "Error: edge pair level already recorded?"
          edge_pair_levels[ep] = j

    #for each edge, build a permutation (list) which says what is required
    #None means there is no requirement
    edge_perms = [ [None for _ in xrange(cover_degree)] for i in xrange(len(self.E))]
    
    for i,p in enumerate(paths):
      for j, (e,d) in enumerate(p):
        j_prev = (j-1)%len(p)
        level = edge_pair_levels[ (i,j) ]
        prev_level = edge_pair_levels[ (i,j_prev) ]
        if d:
          edge_perms[e][j_prev] = j
        else:
          edge_perms[e][j] = j_prev
    
    #now we have partial permutations; we need to fill them in
    #we want to do this in such a way that there are as few loops as possible
    #or each permutations is close to transitive, just to mess things up
    for ei,e in enumerate(self.E):
      source_indices = [i for i in xrange(cover_degree) if edge_perms[ei][i] == None]
      dest_indices = [i for i in xrange(cover_degree) if i not in edge_perms[ei]]
      
      
      

        
        
          
        
      

  
  
  
  def next_edge(self, current_edge, current_direction):
    """returns the next edge (reading along boundary); for directions, 0 means forward, 1 backward"""
    if type(current_direction) == bool:
      current_direction = (0 if current_direction else 1)
    e = self.E[current_edge]
    vert = self.V[(e.dest if current_direction==0 else e.source)]
    to_look_for = (current_edge, (1-current_direction)==0)
    ind = vert.edges.index(to_look_for)
    ind = (ind+1)%len(vert.edges)
    return vert.edges[ind][0], vert.edges[ind][1]  #it's 1- because if it's outgoing (1), then we do it forward (0)
  
  def boundaries(self):
    """returns a list of the boundary words in the fatgraph"""
    edges_read_list = [[False, False] for e in self.E]
    boundaries = []
    while True:
      #find an unread edge direction
      found_edge = False
      for i,e in enumerate(edges_read_list):
        if not e[0]:
          current_edge = i
          current_direction = 0  #0 means forward
          found_edge = True
          break
        elif not e[1]:
          current_edge = i
          current_direction = 1
          found_edge = True
          break
      if not found_edge:
        break
      #now scan forward and read the edges in this boundary
      boundaries.append('')
      while not edges_read_list[current_edge][current_direction]:
        e = self.E[current_edge]
        boundaries[-1] += (e.label_forward if current_direction==0 else e.label_backward)
        edges_read_list[current_edge][current_direction] = True
        current_edge, current_direction = self.next_edge(current_edge, current_direction)
        current_direction = (0 if current_direction else 1)
    return boundaries
  
  def boundaries_with_patterns(self):
    """returns a list of the boundary patterns in a fatgraph"""
    edges_read_list = [[False, False] for e in self.E]
    boundaries = []
    while True:
      #find an unread edge direction
      found_edge = False
      for i,e in enumerate(edges_read_list):
        if not e[0]:
          current_edge = i
          current_direction = 0  #0 means forward
          found_edge = True
          break
        elif not e[1]:
          current_edge = i
          current_direction = 1
          found_edge = True
          break
      if not found_edge:
        break
      #now scan forward and read the edges in this boundary
      boundaries.append([])
      while not edges_read_list[current_edge][current_direction]:
        e = self.E[current_edge]
        outgoing_label = (e.label_forward if current_direction==0 else e.label_backward)
        start_vert_ind = (e.source if current_direction==0 else e.dest)
        all_outgoing_labels = self.outgoing_labels(start_vert_ind)
        to_append = [outgoing_label, [x for x in all_outgoing_labels if x != outgoing_label]]
        boundaries[-1].append(to_append)
        edges_read_list[current_edge][current_direction] = True
        current_edge, current_direction = self.next_edge(current_edge, current_direction)
      #note that we will have one extra outgoing label for every letter 
      #scan through and remove the "outgoing label" which is really the inverse 
      #of the incoming label
      for i in xrange(len(boundaries[-1])):
        ip1 = (i+1)%len(boundaries[-1])
        current_letter = boundaries[-1][i][0]
        olist = boundaries[-1][ip1][-1]
        boundaries[-1][ip1][-1] = [x for x in olist if x != current_letter.swapcase()]
      #put the actual boundary as the first entry
      boundary = ''.join([b[0] for b in boundaries[-1]])
      boundaries[-1] = [boundary] + boundaries[-1]
        
    return boundaries    
  
  def lift(self, G):
    """returns the fatgraph which is self lifted to the finite cover G"""
    new_verts = [Vertex(None, None) for i in xrange(G.degree*len(self.V))]
    new_edges = []
    vertices_over_vert = [range(G.degee*i, G.degree*(i+1)) for i in xrange(len(self.V))]
    for i,v in enumerate(self.V):
      for cover_v_i in vertices_over_vert[i]:
        new_verts[cover_v_i] = copy.deepcopy(v)
    
    #go to every vertex in the cover, and go out every outgoing edge
    #(this will only go over every edge once)
    for i,v in enumerate(self.V):
      for sheet, cover_v_i in enumerate(vertices_over_vert[i]):
        for k, (e, d) in enumerate(v.edges):
          if not d: 
            continue
          source_ind = cover_v_i
          edge_label_f = self.E[e].label_forward
          edge_label_b = self.E[e].label_backward
          if edge_label_f.isupper():
            G_gen_ind = G.base_gen_inds[edge_label_f.swapcase()]
            dest_sheet = G.base_gen_inverse_actions[G_gen_ind][sheet]
          else:
            G_gen_ind = G.base_gen_inds[edge_label_f]
            dest_sheet = G.base_gen_actions[G_gen_ind][sheet]
          dest_base = self.E[e].dest
          dest_ind = vertices_over_vert[dest_base][dest_sheet]
          dest_v_edge_ind = self.V[dest_base].find_edge_ind( (e, not d) )
          new_edge_ind = len(new_edges)
          new_verts[source_ind].edges[k] = (new_edge_ind, d)
          new_verts[dest_ind].edges[dest_v_edge_ind] = (new_edge_ind, not d)
          new_edges.append( Edge(source_ind, dest_ind, edge_label_f, edge_label_b) )
    
    return Fatgraph(new_verts, new_edges)
  
  def ends(self):
    """Returns the list of ends, in cyclic order.  Effectively this just means 
    returning the list of pi_1 generators as words in the free group, 
    but being careful about the fatgraph structure.  This is basically the 
    same function as comb_and_find_gen_words, but it's careful to 
    do a tree search in the right order."""
    self.comb()
    edge_stack = [x for x in self.V[self.basepoint].edges]
    #get all the edge paths from the vertices to the basepoint
    edge_paths = self.all_edge_paths_from_basepoint_to_verts()
    end_list = []
    while len(edge_stack) > 0:
      (ei,d) = edge_stack.pop()
      vi = (self.E[ei].dest if d else self.E[ei].source)
      v = self.V[vi]
      vei = v.find_edge_ind( (ei, not d) )
      nve = len(v.edges)
      i = (vei+1)%nve
      while i != vei:
        (e2i, d2) = v.edges[i]
        if self.E[e2i].towards_basepoint == None:
          #get the path to the current vertex
          path_to_v = edge_paths[vi]
          other_vi = (self.E[e2i].dest if d2 else self.E[e2i].source)
          path_to_other = edge_paths[other_vi]
          path_from_other = [(e3i, not d3) for (e3i, d3) in path_to_other[::-1]]
          end_list.append( path_to_v + [(e2i, not d2)] + path_from_other )
        else:
          edge_stack.append( (e2i, d2) )
        i = (i+1)%nve
    
    #collapse the end list into words
    for i,el in enumerate(end_list):
      label_list = [ (self.E[ei].label_forward if d else self.E[ei].label_backward) for (ei, d) in el]
      end_list[i] = FreeGroupEnd('', ''.join(label_list))
      
    return end_list
  
  
  
def read_file(filename):
  """Read a fatgraph .fg file"""
  f = open(filename, 'r')
  lines = f.read().split('\n')
  #read the vertices, with edge names instead of indices
  current_line = 0
  while lines[current_line][0] == '#':
    current_line += 1
  num_verts = int(lines[current_line].split()[-1])
  current_line += 1
  vert_name_indices = {}
  verts = []
  while len(verts) < num_verts:
    vert_name_indices[lines[current_line].split()[0]] = len(verts)
    new_vertex = Vertex(lines[current_line+1].split(),              \
                        map(int, lines[current_line+2].split()))
    verts.append(new_vertex)
    if lines[current_line+3][:6] == 'bezier':
      current_line += 5
    else:
      current_line += 3
  #and the edges, making the dictionary of edge names -> indices
  num_edges = int(lines[current_line].split()[-1])
  current_line += 1
  edge_name_indices = {}
  edges = []
  while len(edges) < num_edges:
    name, Lf, Lb, source, dest = lines[current_line].split()
    edge_name_indices[name] = len(edges)
    new_edge = Edge(vert_name_indices[source], \
                    vert_name_indices[dest],   \
                    Lf,                        \
                    Lb)                        
    edges.append(new_edge)
    current_line += 1
  
  #refresh the vertex edge names into indices
  for v in verts:
    v.edges = [(edge_name_indices[n], direc) for n,direc in v.edges]
  
  f.close()
  
  return Fatgraph(verts, edges)



def rose_plus_word(g, w):
  """returns a surface of genus g with one boundary component [a,b][c,d]..., one of 
  whose generators is labeled by w"""
  v_edge_list = [ [2*i,2*i+1,2*i,2*i+1] for i in xrange(g)]
  v_edge_list = [x for L in v_edge_list for x in L]
  v_dir_list = g*[True, False, False, True]
  verts = [Vertex(v_edge_list, v_dir_list)]
  edges = [Edge(0,0,alphabet[i], alphabet[i].swapcase()) for i in xrange(2*g-1)]
  for i in xrange(len(w)-1):
    edges.append( Edge(i,i+1,w[i],w[i].swapcase()) )
    verts.append( Vertex([(2*g-1)+i, (2*g-1)+i+1], [False, True]) )
  edges.append( Edge(len(w)-1, 0, w[-1], w[-1].swapcase()) )
  verts[0].edges[-1] = ((2*g-1), True)
  verts[0].edges[-3] = ((2*g-1)+len(w)-1, False)
  return Fatgraph(verts, edges)
  
def fatgraph_from_order(O):
  letters_to_edges = {}
  E = []
  for ell in O:
    if ell.swapcase() in letters_to_edges:
      i,d = letters_to_edges[ell.swapcase()]
      letters_to_edges[ell] = (i, not d)
    else:
      letters_to_edges[ell] = (len(E), ell.islower())
      E.append( Edge(0,0,ell.lower(),ell.upper()) )

  V = [ Vertex([letters_to_edges[ell] for ell in O]) ]
  return Fatgraph(V,E)
  
def free_map_kernel(target_words, verbose=0):
  #build a fatgraph which is just a bunch of loops
  #labeled by the target_words
  V = [ Vertex([]) ]
  E = []
  V[0].is_basepoint = True
  orig_gen_words = []
  #for each word, 
  for w in target_words:
    orig_gen_words.append( [ [], w ] )
    for i in xrange(len(w)):
      this_edge_ind = len(E)
      E.append( Edge( -1, -1, w[i], inverse(w[i]) ) )
      if i == 0:
        prev_vert_ind = 0
        V[0].edges.append( (this_edge_ind, True) )
      else:
        prev_vert_ind = len(V)-1
        V[prev_vert_ind].edges.append( (this_edge_ind, True) )
      if i == len(w)-1:
        next_vert_ind = 0
        V[0].edges.append( (this_edge_ind, False) )
      else:
        next_vert_ind = len(V)
        V.append( Vertex( [ (this_edge_ind, False) ] ) )
      E[-1].source = prev_vert_ind
      E[-1].dest = next_vert_ind
      orig_gen_words[-1][0].append( (this_edge_ind, True) )
    E[-1].gen_word = [(len(orig_gen_words)-1, 1)]
  G = Fatgraph(V,E)
  G.orig_gen_words = orig_gen_words
  G.basepoint = 0
#  return G
  G = G.fold(verbose=(verbose>0))
  #the kernel words will be in G.kernel_words
  #the source gens will be x,y,z,w if there's <= 4, or x1, x2, ... if more
  if len(target_words) <= 4:
    gens = ['x','y','z','w'][:len(target_words)]
  else:
    gens = ['x' + str(i) for i in xrange(len(target_words))]
  ker_words = []
  for kw in G.kernel_words:
    ker_words.append('')
    for (g,s) in kw:
      ker_words[-1] += (gens[g] if s>0 else inverse(gens[g]))
  return ker_words






















