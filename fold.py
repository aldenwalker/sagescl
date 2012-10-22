#!/usr/bin/python


class Vertex:
  def __init__(self, p, n):
    self.edges = []
    self.next = n
    self.prev = p
    if n != None:
      n.prev = self
    if p != None:
      p.next = self

  def kill_me(self):
    if self.next == None:
      if self.prev == None:
        return
      self.prev.next = None
      return
    if self.prev == None:
      self.next.prev = None
      return
    self.next.prev = self.prev
    self.prev.next = self.next
    return
      

class Edge:
  def __init__(self, L, v1, v2, p, n):
    self.next = n
    self.prev = p
    if n != None:
      n.prev = self
    if p != None:
      p.next = self
    self.label = L
    self.start = v1
    self.end = v2
    v1.edges.append( (self, 'out') )
    v2.edges.append( (self, 'in') )
    
  def kill_me(self):
    if self.next == None:
      if self.prev == None:
        return
      self.prev.next = None
      return
    if self.prev == None:
      self.next.prev = None
      return
    self.next.prev = self.prev
    self.prev.next = self.next
    return

class Graph:
  def __init__(self, words=None):
    self.first_edge = None
    self.first_vertex = None
    self.origin_vertex = None
    if words == None:
      return
    self.first_vertex = Vertex(None, None)
    self.origin_vertex = self.first_vertex
    for w in words:
      if len(w) > 1:
        self.first_vertex = Vertex(None, self.first_vertex)
      self.first_edge = Edge( w[0], self.origin_vertex, self.first_vertex, None, self.first_edge)
      for i in xrange(1, len(w)-1):
        self.first_vertex = Vertex(None, self.first_vertex)
        self.first_edge = Edge( w[i], self.first_vertex.next, self.first_vertex, None, self.first_edge) 
      if len(w) > 1:
        self.first_edge = Edge( w[-1], self.first_vertex, self.origin_vertex, None, self.first_edge) 
  
  def add_words(self, words) :
    for w in words:
      if len(w) > 1:
        self.first_vertex = Vertex(None, self.first_vertex)
      self.first_edge = Edge( w[0], self.origin_vertex, self.first_vertex, None, self.first_edge)
      for i in xrange(1, len(w)-1):
        self.first_vertex = Vertex(None, self.first_vertex)
        self.first_edge = Edge( w[i], self.first_vertex.next, self.first_vertex, None, self.first_edge)
      if len(w) > 1:
        self.first_edge = Edge( w[-1], self.first_vertex, self.origin_vertex, None, self.first_edge) 

  
    
  def __str__(self):
    ans = ''
    ans += "Vertices:\n"
    v = self.first_vertex
    while v != None:
      ans += str(v) + ', ' + str([e[0] for e in v.edges]) + '\n'
      v = v.next
    ans += 'Edges:\n'
    e = self.first_edge
    while e != None:
      ans += str(e) + ',' + e.label + ' ' + str(e.start) + ' -> ' + str(e.end) +'\n'
      e = e.next
    return ans
  
  #def cpy(self):
  #  new_graph = Graph()
  #  new_graph.vertices = [v for v in self.vertices]
  #  new_graph.edges = [e for e in self.edges]
  #  for v in new_graph.vertices:
  #    v.edges = [ ( new_graph.edges[ e[0].ind ], e[1]) for e in v.edges]
  #  for e in new_graph.edges:
  #    e.start = new_graph.vertices[e.start.ind]
  #    e.end = new_graph.vertices[e.end.ind]
  #  return new_graph
    
  def unfolded_vertex(self):
    v = self.first_vertex
    while v != None:
      for j in xrange(len(v.edges)):
        for k in xrange(j+1, len(v.edges)):
          if v.edges[j][0].label == v.edges[k][0].label and \
             v.edges[j][1] == v.edges[k][1]:   # if the labels are the same and same direction
            return (v,j,k)
          if v.edges[j][0].label == v.edges[k][0].label.swapcase() and \
             v.edges[j][1] != v.edges[k][1]:   # labels are inverse and different directions
            return (v,j,k)
      v = v.next
    return None
      
  def fold(self):
    num_folds = 0
    while True:
      #print str(self)
      i = self.unfolded_vertex()
      if i == None:
        return num_folds
      v, j, k = i
      
      #print "Found unfolded edges ", j, " and ", k, " in vertex ", v
      
      if v.edges[j][1] == 'in':
        vertex_1 = v.edges[j][0].start 
        dir_1 = 'out'
      else:
        vertex_1 = v.edges[j][0].end
        dir_1 = 'in'
        
      if v.edges[k][1] == 'in':
        vertex_2 = v.edges[k][0].start 
        dir_1 = 'out'
      else:
        vertex_2 = v.edges[k][0].end
        dir_1 = 'in'
      
      #if the vertices are equal, then we remove the edge
      if vertex_1 == vertex_2:
        #print "Removing edge ", self.vertices[i].edges[k][0].ind
        e_to_kill = v.edges[k][0]
        vertex_1.edges.remove( (e_to_kill, dir_1) )
        v.edges.remove( (e_to_kill, ('in' if dir_1 == 'out' else 'out')) )
        #print "Deleting edge ", e_to_kill.ind
        e_to_kill.kill_me()
        #del e_to_kill
        num_folds += 1
      #if the vertices aren't equal, then make them equal
      #note don't remove the edge yet, just for simplicity
      else:
        for m in xrange(len(vertex_2.edges)):
          if vertex_2.edges[m][0].start == vertex_2:
            vertex_2.edges[m][0].start = vertex_1
          if vertex_2.edges[m][0].end == vertex_2:
            vertex_2.edges[m][0].end = vertex_1
        vertex_1.edges.extend(vertex_2.edges)
        vertex_2.kill_me()
        #del vertex_2
        num_folds -= 1  # note we *will* kill a loop because of the vertex we just killed, but that doesn't count























