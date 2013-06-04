"""This module implements a fatgraph class"""

import copy

class Edge:
  def __init__(self, v0, v1, Lf, Lb):
    self.source = v0
    self.dest = v1
    self.label_forward = Lf
    self.label_backward = Lb
    self.carries_folded_edges = None
    self.carried_by_edge = None
  
  def __str__(self):
    s = '(' + str(self.source) + '->' + str(self.dest) + ', ' + self.label_forward + ', ' + self.label_backward + ')'
    if self.carries_folded_edges != None:
      s += '; folded: ' + str(self.carries_folded_edges)
    if hasattr(self, 'dead'):
      s += "*dead*"
    if self.carried_by_edge != None:
      s += ' carried by ' + str(self.carried_by_edge)
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
    self.carries_folded_verts = None
    self.carried_by_vert = None
    
  def __str__(self):
    s =  str(self.edges) 
    if self.carries_folded_verts != None:
      s += '; folded: ' + str(self.carries_folded_verts)
    if hasattr(self, 'dead'):
      s += "*dead*"
    if self.carried_by_vert != None:
      s += ' carried by ' + str(self.carried_by_vert)
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
      
    return ans
  
  def outgoing_labels(self, ind):
    edges = self.V[ind].edges
    ans = []
    for i,direc in edges:
      if direc:  #it's outgoing
        ans.append( self.E[i].label_forward )
      else:
        ans.append( self.E[i].label_backward )
    return ans
  
  def cleanup(self):
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
  
  
  def unfolded_edge_pair(self, v_ind):
    """returns a pair of edge indices which have the same outgoing label,
       or None if no such pair exists"""
    outgoing_labels = {}
    for i in xrange(len(self.V[v_ind].edges)):
      edge_ind, edge_dir = self.V[v_ind].edges[i]
      ol = (self.E[edge_ind].label_forward if edge_dir else self.E[edge_ind].label_backward)
      if ol in outgoing_labels:
        return (outgoing_labels[ol], i)
      outgoing_labels[ol] = i
    return None
  
  def fold(self): 
    """returns the folded version of the fatgraph, with the folded structure"""
    #initialize the new folded fatgraph
    new_F = Fatgraph(self.V, self.E)
    new_F.unfolded_V = copy.deepcopy(self.V)
    new_F.unfolded_E = copy.deepcopy(self.E)
    for i in xrange(len(new_F.V)):
      new_F.V[i].carries_folded_verts = [i]
    for i in xrange(len(new_F.E)):
      new_F.E[i].carries_folded_edges = [i]
    
    while True:
      print "Current fatgraph: "
      print new_F
      
      #find an unfolded vertex
      unfolded_vert = None
      unfolded_e_p = None #which particular pair of edge indices in the vert are duplicates
      for i in xrange(len(new_F.V)):
        unfolded_e_p = new_F.unfolded_edge_pair(i)
        if unfolded_e_p != None:
          unfolded_vert = i
          break
      if unfolded_vert == None:
        break
      
      print "Found unfolded vertex ", unfolded_vert, " with edges ", unfolded_e_p
      
      #fold the edges together
      i1, i2 = unfolded_e_p
      v = new_F.V[unfolded_vert]
      e1, e2 = new_F.E[v.edges[i1][0]], new_F.E[v.edges[i2][0]]
      e1.carries_folded_edges.extend(e2.carries_folded_edges)
      e2.dead = True
      
      print "This is edges with main indices ", v.edges[i1], ' and ', v.edges[i2]
      
      #fold the vertices together    
      other_vert_1 = (e1.dest if v.edges[i1][1] else e1.source)
      other_vert_2 = (e2.dest if v.edges[i2][1] else e2.source)
      
      print "The vertices to fold together are ", other_vert_1, ' and ', other_vert_2
      
      other_vert_2_ind = new_F.V[other_vert_2].find_edge_ind(v.edges[i2][0], not v.edges[i2][1])
      if unfolded_vert == other_vert_2:
        #we need to erase both v.edges[i2] and v.edges[other_vert_2_ind]
        del v.edges[max(i2, other_vert_2_ind)]
        del v.edges[min(i2, other_vert_2_ind)]
      else:
        del v.edges[i2]
        del new_F.V[other_vert_2].edges[other_vert_2_ind]
        
      if other_vert_1 != other_vert_2:
        for i in xrange(len(new_F.V[other_vert_2].edges)):
          ei = new_F.V[other_vert_2].edges[i]
          if ei[1]:
            new_F.E[ei[0]].source = other_vert_1
          else:
            new_F.E[ei[0]].dest = other_vert_1
        new_F.V[other_vert_1].edges.extend(new_F.V[other_vert_2].edges)
        new_F.V[other_vert_1].carries_folded_verts.extend(new_F.V[other_vert_2].carries_folded_verts)
        new_F.V[other_vert_2].dead = True
    
    #the graph is folded now, but we need to clean it up
    new_F.cleanup()
    return new_F
    
  
  def non_injective_pieces(self):
    """for a folded fatgraph, returns the rectangles and polygons 
    which can be assembled to give any loop in the kernel"""
    
    rectangles = []
    polygons = []
    
    #go through all the edges, and look at the original edges they carry
    #any pair of original edges they carry is a legitimate rectangle
    #a rectangle is a pair ((e1, bool), (e2,bool)), where the bools
    #record whether the edge is forward or backward
    
    for e in self.E:
      for i in xrange(len(e.carries_folded_edges)):
        e1 = e.carries_folded_edges[i]
        dir1 = (self.unfolded_E[e1].forward_label == e.forward_label)
        for j in xrange(i+1, len(e.carries_folded_edges)):
          e2 = e.carries_folded_edges[j]
          dir2 = (self.unfolded_E[e2].forward_label == e.forward_label)
          
  
  
  
  def next_edge(self, current_edge, current_direction):
    """returns the next edge (reading along boundary); for directions, 0 means forward, 1 backward"""
    if type(current_direction) == bool:
      current_direction = (0 if current_direction else 1)
    e = self.E[current_edge]
    vert = self.V[(e.dest if current_direction==0 else e.source)]
    to_look_for = (current_edge, current_direction==0)
    ind = vert.edges.index(to_look_for)
    ind = (ind+1)%len(vert.edges)
    return vert.edges[ind][0], not vert.edges[ind][1]  #it's 1- because if it's outgoing (1), then we do it forward (0)
  
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