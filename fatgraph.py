"""This module implements a fatgraph class"""


from word import *
from sage.all import *

import copy


def approx_rat(x, tol=0.0000001):
  """returns a rational approximation which is closer than tol"""
  c = continued_fraction_list(x, partial_convergents=True, nterms=10)
  for r in c[1]:
    if abs(Integer(r[0])/r[1] - x) < tol:
      return Integer(r[0])/r[1]
  return None
    



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
      else:
        self.edges = []
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
          return (i, (i+1)%lve)
        ov_e1_ind = self.V[ovi1].edges.index( (e1, not d1) )
        ov_e2_ind = self.V[ovi1].edges.index( (e2, not d2) )
        if ov_e1_ind == (ov_e2_ind+1)%len(self.V[ovi1].edges):
          return (i, (i+1)%lve)
    return None
  
  def fold(self, fatgraph_fold=False, verbose=False): 
    """returns the folded version of the fatgraph, with the folded structure; 
    if fatgraph_fold=True, only do fatgraph (surface) folds"""
    #initialize the new folded fatgraph
    new_F = Fatgraph(copy.deepcopy(self.V), copy.deepcopy(self.E))
    new_F.unfolded_V = copy.deepcopy(self.V)
    new_F.unfolded_E = copy.deepcopy(self.E)
    for i in xrange(len(new_F.V)):
      new_F.V[i].carries_folded_verts = [i]
    for i in xrange(len(new_F.E)):
      new_F.E[i].carries_folded_edges = [i]
    
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
          unfolded_e_p = new_F.unfolded_edge_pair(i)
        if unfolded_e_p != None:
          unfolded_vert = i
          break
      if unfolded_vert == None:
        break
      
      if verbose:
        print "Found unfolded vertex ", unfolded_vert, " with edges ", unfolded_e_p
      
      #fold the edges together
      i1, i2 = unfolded_e_p
      v = new_F.V[unfolded_vert]
      ei1, d1 = v.edges[i1]
      ei2, d2 = v.edges[i2]
      e1, e2 = new_F.E[ei1], new_F.E[ei2]
      e1.carries_folded_edges.extend(e2.carries_folded_edges)
      e2.dead = True
      
      if verbose:
        print "This is edges with main indices ", v.edges[i1], ' and ', v.edges[i2]
      
      #get the other vertices 
      ovi1 = (e1.dest if v.edges[i1][1] else e1.source)
      ovi2 = (e2.dest if v.edges[i2][1] else e2.source)
      ov1, ov2 = new_F.V[ovi1], new_F.V[ovi2]
      
      if verbose:
        print "The vertices to fold together are ", ovi1, ' and ', ovi2
      
      #remove edge 2 altogether
      #note we remove the origin first, *then* the destination, and we 
      #remember the destination index.  This ensure we know where to glue 
      #the edges from vertex 1
      del v.edges[i2]
      ov2_e_ind = ov2.edges.index( (ei2, not d2) )
      
      if ovi1 != ovi2:
        #get a list of the edges that are in vertex 2, and in the 
        #correct order!
        edges_from_vert_2 = ov2.edges[ov2_e_ind+1:] \
                          + ov2.edges[:ov2_e_ind]
        
        #for all these edges, make sure they point to the right place
        for (et, dt) in edges_from_vert_2:
          if dt:
            new_F.E[et].source = ovi1
          else:
            new_F.E[et].dest = ovi1
        
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
        del ov2.edges[ov2_e_ind]
    
    #the graph is folded now, but we need to clean it up
    new_F.cleanup()
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
  
  
  def non_injective_pieces(self):
    
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
  
  def kernel_element(self, required_vert=None, verbose=0):
    """return an element in the kernel"""
    #we do a linear programming problem
    R, GE, T = self.non_injective_pieces()
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
    
    #for each edge, add the constraint
    #notice that every triangle and rectangle cannot contain an edge more than once
    #so it contributes only -1, 0,or 1
    for e in GE:
      lf = 0
      eb = (e[1], e[0])
      for i,r in enumerate(R):
        rbound = self.rectangle_boundary(r)
        if e in rbound:
          lf += x[i]
        elif eb in rbound:
          lf += (-1)*x[i]
      for i,t in enumerate(T):
        if e in t:
          lf += x[num_rects + i]
        elif eb in t:
          lf += (-1)*v[num_rects + i]
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
      print "Done creating LP."
      sys.stdout.flush()
      if verbose > 1:
        p.show()
     
    try:
      obj_L1 = p.solve()
    except :
      print "LP Error (probably no kernel)"
      return None
    x_solution = p.get_values(x)
    
    
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
    have_visted_piece = [False for i in xrange(num_rect_pieces)]
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
          next_ind_in_piece = next_p.index( outgoing_edge_inv )
        
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
    """returns a fatgraph which covers self and which contains a lift of path 
    which is embedded, also returns the lifted path"""
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
  
  




















