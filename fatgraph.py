"""This module implements a fatgraph class"""

class Edge:
  def __init__(self, v0, v1, Lf, Lb):
    self.source = v0
    self.dest = v1
    self.label_forward = Lf
    self.label_backward = Lb
    self.carries_folded_edges = None
  
  def __str__(self):
    return '(' + str(self.source) + '->' + str(self.dest) + ', ' + self.label_forward + ', ' + self.label_backward + ')'
    
  def __repr__(self):
    return 'Edge(' + ','.join(map(str, [self.source, self.dest, self.label_forward, self.label_backward])) + ')'
  
class 
  
class Vertex:
  """a vertex class with ordered edges; note "true" means the edge is leaving"""
  def __init__(self, edge_list, direction_list):
    self.edges = [(e, direction_list[i]==1) for i,e in enumerate(edge_list)]
    self.carries_folded_verts = None
    
  def __str__(self):
    return str(self.edges)
  
  def __repr__(self):
    return str(self)
  
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
    self.folded_V = None
    self.folded_E = None
  
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
    
    
    
    #find an unfolded vertex
    unfolded_vert = None
    unfolded_e_p = None #which particular pair of edge indices in the vert are duplicates
    for i in xrange(len(self.V)):
      unfolded_e_p = self.unfolded_edge_pair(i)
      if unfolded_e_p != None:
        unfolded_vert = i
        break
    
  
  def next_edge(self, current_edge, current_direction):
    """returns the next edge (reading along boundary); for directions, 0 means forward, 1 backward"""
    e = self.E[current_edge]
    vert = self.V[(e.dest if current_direction==0 else e.source)]
    to_look_for = (current_edge, current_direction)
    ind = vert.edges.index(to_look_for)
    ind = (ind+1)%len(vert.edges)
    return vert.edges[ind][0], 1-vert.edges[ind][1]  #it's 1- because if it's outgoing (1), then we do it forward (0)
  
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