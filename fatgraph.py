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
  
  def fold(self, verbose=False): 
    """returns the folded version of the fatgraph, with the folded structure"""
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
      e1, e2 = new_F.E[v.edges[i1][0]], new_F.E[v.edges[i2][0]]
      e1.carries_folded_edges.extend(e2.carries_folded_edges)
      e2.dead = True
      
      if verbose:
        print "This is edges with main indices ", v.edges[i1], ' and ', v.edges[i2]
      
      #fold the vertices together    
      other_vert_1 = (e1.dest if v.edges[i1][1] else e1.source)
      other_vert_2 = (e2.dest if v.edges[i2][1] else e2.source)
      
      if verbose:
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
    
    #first let's get a list of the incoming and outgoing letters for 
    #every vertex
    incoming_letters = len(self.unfolded_V)*[[]]
    outgoing_letters = len(self.unfolded_V)*[[]]
    for i,v in enumerate(self.unfolded_V):
      outgoing_letters[i] = [(self.unfolded_E[ei].label_forward if edir else self.unfolded_E[ei].label_backward) for ei, edir in v.edges]
      outgoing_letters[i] = list(set(outgoing_letters[i]))
    incoming_letters = map(inverse, outgoing_letters)
    
    rectangles = []
    
    #go through all the edges, and look at the original edges they carry
    #any pair of original edges they carry is a legitimate rectangle
    #a rectangle is a pair ((e1, bool, w1), (e2,bool, w2)), where the bools
    #record whether the edge is forward or backward, and the w's record 
    #the three-letter word centered on the rectangle
        
    for e in self.E:
      for i in xrange(len(e.carries_folded_edges)):
        e1 = e.carries_folded_edges[i]
        dir1 = (self.unfolded_E[e1].label_forward == e.label_forward)
        e1letter = e.label_forward
        e1_ivert = (self.unfolded_E[e1].source if dir1 else self.unfolded_E[e1].dest)
        e1_dvert = (self.unfolded_E[e1].dest if dir1 else self.unfolded_E[e1].source)
        preceeding_letter_ops1 = [ell for ell in incoming_letters[e1_ivert] if ell != inverse(e1letter)]
        succeeding_letter_ops1 = [ell for ell in outgoing_letters[e1_dvert] if ell != inverse(e1letter)]
        for j in xrange(i+1, len(e.carries_folded_edges)):
          e2 = e.carries_folded_edges[j]
          dir2 = (self.unfolded_E[e2].label_forward != e.label_forward)
          e2letter = e.label_backward
          e2_ivert = (self.unfolded_E[e2].source if dir2 else self.unfolded_E[e2].dest)
          e2_dvert = (self.unfolded_E[e2].dest if dir2 else self.unfolded_E[e2].source)
          preceeding_letter_ops2 = [ell for ell in incoming_letters[e2_ivert] if ell != inverse(e2letter)]
          succeeding_letter_ops2 = [ell for ell in outgoing_letters[e2_dvert] if ell != inverse(e2letter)]
          #print (j,e2,dir2,e2letter,e2_ivert, e2_dvert)
          if e1_ivert == e2_dvert and e1_dvert == e2_ivert:
            #this is a totally dummy rectangle; I guess add it twice
            rectangles.append( ((e1, dir1, e1letter), (e2, dir2, e2letter)) )
            #rectangles.append( ((e2, not dir2, inverse(e2letter)), (e1, not dir1, inverse(e1letter))) )
          elif e1_ivert == e2_dvert:
            #the first edge is a dummy edge; we only need to iterate over possible other letters
            for s1 in succeeding_letter_ops1:
              for p2 in preceeding_letter_ops2:
                rectangles.append( ((e1, dir1, e1letter+s1), (e2, dir2, p2+e2letter)) )
                rectangles.append( ((e2, not dir2, inverse(p2+e2letter)), (e1, not dir1, inverse(e1letter+s1))) )
          elif e1_dvert == e2_ivert:
            #the second edge is a summy edge; we only need to do the other letters
            for p1 in preceeding_letter_ops1:
              for s2 in succeeding_letter_ops2:
                rectangles.append( ((e1, dir1, p1+e1letter), (e2, dir2, e2letter+s2)) )
                rectangles.append( ((e2, not dir2, inverse(e2letter+s2)), (e1, not dir1, inverse(p1+e1letter))) )
          else:
            #both edges will be nontrivial; we need to look at all of them
            for p1 in preceeding_letter_ops1:
              for s1 in succeeding_letter_ops1:
                for p2 in preceeding_letter_ops2:
                  for s2 in succeeding_letter_ops2:
                    rectangles.append( ((e1, dir1, p1+e1letter+s1), (e2, dir2, p2+e2letter+s2)) )
                    rectangles.append( ((e2, not dir2, inverse(p2+e2letter+s2)), (e1, not dir1, inverse(p1+e1letter+s1))) )
    
    #an edge is ( (v1, w1), (v2, w2) ), where the vi are unfolded vertices, 
    #and the wi are pairs of letters which are incoming, outgoing
    #the inverse of edge ((v1,w1),(v2,w2)) is ((v2,w2),(v1,w1))
    #note, if v1 == v2, then we can cut the tree to reduce its size, so 
    #we may assume that v1 != v2
    #we will assume that v1 < v2 for uniqueness
    #also note, v1 and v2 have to be carried by the same folded vertex
    gluing_edges = []
    ws = [[w11+w12 for w11 in incoming_letters[i] for w12 in outgoing_letters[i] if w12 != w11.swapcase()] for i in xrange(len(self.unfolded_V))]
    for i,v in enumerate(self.V):
      for j in v.carries_folded_verts:
        for k in v.carries_folded_verts:
          if k <= j:
            continue
          for w1 in ws[j]:
            for w2 in ws[k]:
              gluing_edges.append( ((j,w1),(k,w2)) )
    
    #this is an argument not to use C++:
    gluing_edge_index = dict([ (gluing_edges[i], i) for i in xrange(len(gluing_edges))])
    
    #now we must build the triangles
    #this is simply all triangles (triple of distinct vertices-plus-words)
    #a triangle is a triple ((v1,w1), (v2, w2), (v2,w3))
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
            for w1 in ws[i]:
              for w2 in ws[j]:
                for w3 in ws[k]:
                  triangles.append( ((i,w1),(j,w2),(k,w3)) )
    
    return rectangles, gluing_edges, triangles
    
  def rectangle_boundary(self, r):
    """returns the boundary of a rectangle, as a tuple of two edges"""
    ((e1, dir1, w1), (e2, dir2, w2)) = r
    v11 = (self.unfolded_E[e1].source if dir1 else self.unfolded_E[e1].dest)
    v12 = (self.unfolded_E[e1].dest if dir1 else self.unfolded_E[e1].source)
    v21 = (self.unfolded_E[e2].source if dir2 else self.unfolded_E[e2].dest)
    v22 = (self.unfolded_E[e2].dest if dir2 else self.unfolded_E[e2].source)

    edge1 = ( ((v22,w2[-2:]),(v11,w1[:2])) if v11 != v22 else (None, v11))
    edge2 = ( ((v12,w1[-2:]),(v21,w2[:2])) if v21 != v12 else (None, v21))
    return (edge1, edge2)
  
  def rectangle_verts(self, r):
    """returns a list (potentially with dupes) of the vertices in a rectangle"""
    ((e1, dir1, w1), (e2, dir2, w2)) = r
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
      edge1, edge2 = self.rectangle_boundary(r)
      if edge1[0] == None:
        coef += 1
      if edge2[0] == None:
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
     
    obj_L1 = p.solve()
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
        if outgoing_edge[0] == None:
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

    
    
  
  
  def kernel_elements_old(self):
    """return an element in the kernel of the map to the free group"""
    R = self.non_injective_rectangles()
    
    #the stack keeps all the current tree pieces, all of which should have one 
    #boundary edge.  this starts as just the rectangles which pinch off
    
    S = []
    
    for i in xrange(len(R)):
      (e1, dir1), (e2, dir2) = R[i]
      E1, E2 = self.unfolded_E[e1], self.unfolded_E[e2]
      iv1, dv1 = ((E1.source, E1.dest) if dir1 else (E1.dest, E1.source))
      iv2, dv2 = ((E2.source, E2.dest) if dir2 else (E2.dest, E2.source))
      if self.unfolded_V[iv1].carried_by_vert != self.unfolded_V[dv2].carried_by_vert:
        print "Error edges don't make sense"
      elif self.unfolded_V[dv1].carried_by_vert != self.unfolded_V[iv2].carried_by_vert:
        print "Error edges don't make sense"
      if iv1 == dv2:
        #start a fatgraph
        e = Edge(0,1,(e1, dir1), (e2, dir2))
        v0 = Vertex([(0, True)])
        v1 = Vertex([(0, False)])
        F = Fatgraph([v01, v1], [e])
        F.open_position = (1, 0) #1st vertex, 0th position
        F.open_edge = ( (dv1, (E1.label_forwards if dir1 else E1.label_backward)), \
                        (iv2, (E2.label_forwards if dir2 else E2.label_backward)) )
        F.open_vertex = self.unfolded_V[dv1].carried_by_vert
        F.open_vertex_arrival_time = 0
        S.append(F)
      if dv1 == iv2:
        e = Edge(0,1,(e2, dir2), (e1, dir1))
        v0 = Vertex([(0, True)])
        v1 = Vertex([(0, False)])
        F = Fatgraph([v01, v1], [e])
        F.open_position = (1,0)
        F.open_edge = ( (dv2, (E2.label_forwards if dir2 else E2.label_backward)), \
                        (iv1, (E1.label_forwards if dir1 else E1.label_backward)) )
        F.open_vertex = self.unfolded_V[dv2].carried_by_vert
        F.open_vertex_arrival_time = 0
        S.append(F)
        
    # we've initialized the set of trees 
    # at every step, there are two stages.
    #
    # First, we see if any open edges can be simply extended
    # 
    # next, we collect the trees that are at each vertex
    # into strings and loops -- a loop gives a complete tree which is finished
    # a loop gives a new tree with an open edge at the same vertex
    # a string or loop must have some tree with arrival time now 
    # (otherwise, we'd generate the same trees multiple times)
    while True:
      # collect trees into strings and loops
      trees_at_vert = {}
      for i in xrange(len(S)):
        trees_at_vert[S[i].open_vertex] = trees_at_vert.get(S[i].open_vertex, []) + [i]
      for ind in trees_at_vert:
        tav = trees_at_vert[ind]
        
  
  
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
  





















