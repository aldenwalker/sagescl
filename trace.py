#!/usr/bin/python


import subprocess
import copy
import os
#import mpmath


FNULL = open(os.devnull, 'w')
#CWD = os.getcwd()
#HOME = os.environ["HOME"]
QPATH = '/home/akwalker/Documents/software/qhull-2012.1/bin/'

def multiply_words(w1, w2=None):
  if type(w1) == str:
    i = 0
    w1L = len(w1)
    w2L = len(w2)
    while i < w1L and i < w2L and w1[w1L-i-1] == w2[i].swapcase():
      i += 1
    return w1[:w1L-i] + w2[i:]
  elif type(w1) == list:
    return reduce(multiply_words, w1)
    
def inverse(w):
  return w[::-1].swapcase()

def commutator(x,y):
  return multiply_words([x,y,inverse(x),inverse(y)])


def convex_hull(L):
  inp = '{0}\n{1}\n'.format(len(L[0]), len(L))
  inp = inp + '\n'.join([ '{0} {1} {2}'.format(*l) for l in L])

  #print inp
  
  #print QPATH+'qhull'
  
  qh = subprocess.Popen([QPATH+'qhull', 'p'],  \
                         stdin=subprocess.PIPE,                    \
                         stdout=subprocess.PIPE,                  \
                         stderr=FNULL)
  qhout = qh.communicate(inp)[0]
  
  #print qhout
  #print ''
  
  s = qhout.split('\n')
  numVerts = int(s[1])
  verts = s[2:(2+numVerts)]
  verts = [map(int, x.split()) for x in verts]
  return verts


def convex_hull_facets(L):
  inp = '{0}\n{1}\n'.format(len(L[0]), len(L))
  inp = inp + '\n'.join([ '{0} {1} {2}'.format(*l) for l in L])

  #print inp

  qh = subprocess.Popen([QPATH+'qhull', 'n'],  \
                         stdin=subprocess.PIPE,                    \
                         stdout=subprocess.PIPE,                  \
                         stderr=FNULL)
  qhout = qh.communicate(inp)[0]
  print qhout
  lines = qhout.split('\n')
  #print lines
  normals = lines[2:]
  normals = [x.split() for x in normals if x!='']
  normals = [map(float, x) for x in normals]
  print normals
  
  #ok now go through the vertices and group them into facets
  facets = []
  for n in normals:
    facets.append([])
    for vert in L:
      d = vert[0]*n[0] + vert[1]*n[1] + vert[2]*n[2]
      if abs(d-(-n[3])) < 1e-6:
        facets[-1].append(vert)
  
  #now we have the facets; find the edges
  numFacets = len(facets)
  edges = []
  for i in xrange(numFacets):
    for j in xrange(i+1, numFacets):
      edges.append([])
      for vert in L:
        d1 = vert[0]*normals[i][0] + vert[1]*normals[i][1] + vert[2]*normals[i][2]
        d2 = vert[0]*normals[j][0] + vert[1]*normals[j][1] + vert[2]*normals[j][2]
        if abs(d1-(-normals[i][3])) < 1e-6 and abs(d2-(-normals[j][3])) < 1e-6:
          edges[-1].append(vert)
  
  edges = [x for x in edges if x != []]
  
  return (facets, edges)






#a polynomial class for multiple variables 
class poly:
  def __init__(self, *args):
    self.vars = sorted(args)
    self.degree = -1
    self.data = {}
    
  def addVars(self, varsToAdd):
    newVars = sorted(self.vars + varsToAdd)
    newLen = len(newVars)
    oldLen = len(self.vars)
    newIndices = [newVars.index(x) for x in self.vars]
    newData = {}
    for term in self.data:
      newTuple = [0 for x in xrange(newLen)]
      for i in xrange(oldLen):
        newTuple[newIndices[i]] = term[i]
      newData[tuple(newTuple)] = self.data[term]
    self.data = newData
    self.vars = newVars
    
  def delVars(self, varsToDel):
    oldLen = len(self.vars)
    if len(varsToDel) == oldLen:
      print "Can't delete all vars"
    newLen = oldLen - len(varsToDel)
    newData = {}
    goodIndices = [x for x in xrange(oldLen) if self.vars[x] not in varsToDel]
        
    for y in varsToDel:
      self.vars.remove(y)
      
    for term in self.data:
      newTuple = [term[i] for i in xrange(oldLen) if i in goodIndices]
      newData[tuple(newTuple)] = self.data[term]
      
    self.data = newData
      
    
  def cleanup(self):
    for k in copy.deepcopy(self.data.keys()):
      if self.data[k] == 0 and not all([d == 0 for d in k]):
        del self.data[k]
    varsToDel = []
    for i in xrange(len(self.vars)):
      if all( [x[i] == 0 for x in self.data] ):
        varsToDel.append(self.vars[i])
    if len(varsToDel) == len(self.vars):
      varsToDel = varsToDel[:-1]
    self.delVars(varsToDel)
    
    
  def __add__(self, other):
    if self.vars != other.vars:
      poly1 = copy.deepcopy(self)
      poly1.addVars([x for x in other.vars if x not in self.vars])
      poly2 = copy.deepcopy(other)
      poly2.addVars([x for x in self.vars if x not in other.vars])
    else:
      poly1 = self
      poly2 = other
    ans = poly(*poly1.vars)
    for term in poly1.data:
      ans.data[term] = ans.data.get(term, 0) + poly1.data[term]
    for term in poly2.data:
      ans.data[term] = ans.data.get(term, 0) + poly2.data[term]
    ans.cleanup()
    ans.degree = max([sum(x) for x in ans.data])
    return ans
  
  def __neg__(self):
    a = copy.deepcopy(self)
    for x in a.data:
      a.data[x] *= -1
    return a
    
  def __sub__(self, other):
    return self + (-other)
    
  def __mul__(self, other):
    if self.vars != other.vars:
      poly1 = copy.deepcopy(self)
      poly1.addVars([x for x in other.vars if x not in self.vars])
      poly2 = copy.deepcopy(other)
      poly2.addVars([x for x in self.vars if x not in other.vars])
    else:
      poly1 = self
      poly2 = other
    ans = poly(*poly1.vars)
    ans.degree = poly1.degree + poly2.degree
    nv = len(poly1.vars) 
    for term in poly1.data:
      for term2 in poly2.data:
        termInAnswer = tuple([term[i]+term2[i] for i in xrange(nv)])
        ans.data[termInAnswer] = ans.data.get(termInAnswer, 0) + (poly1.data[term] * poly2.data[term2])
    ans.cleanup()
    return ans
    
  #this is not as fast as it could be
  def __pow__(self, ex):
    if ex == 0:
      p = poly('x')
      p.degree = 0
      p.data = {tuple([0]):1}
      return p
    elif ex == 1:
      return copy.deepcopy(self)
    elif ex%2 == 0:
      p = self ** (ex/2)
      return p*p
    else:
      return self * (self ** (ex-1))
    
  def __repr__(self):
    return "Polynomial in " + ', '.join(self.vars) + " of degree " + str(self.degree) + ':\n' + str(self)
    
  def tupleToStr(self, tup):
    var_string = [ '(' + self.vars[i] + ')' + ('^' + str(tup[i]) if tup[i]>1 else '') for i in xrange(len(tup))]
    return '*'.join(var_string)
  
  def __str__(self):
    if len(self.data) == 0:
      return ''
    sk = sorted(self.data.keys())
    ans = ''
    
    if len(sk) > 1 and all([w==0 for w in sk[0]]) and self.data[sk[0]]==0:
      del sk[0]
    
    if all([w==0 for w in sk[0]]):
        ans += str(self.data[sk[0]])
    else:
      if self.data[sk[0]] < 0:
        ans += ((str(self.data[sk[0]])+'*') if self.data[sk[0]] != -1 else '-') + self.tupleToStr(sk[0])
      else:
        ans += ((str(self.data[sk[0]])+'*') if self.data[sk[0]] != 1 else '') + self.tupleToStr(sk[0])
    for term in sk[1:]:
      if self.data[term] < 0:
        ans +=  ' - ' + ((str(-self.data[term])+'*') if self.data[term] != -1 else '') + self.tupleToStr(term)
      else:
        ans += ' + ' + ((str(self.data[term])+'*') if self.data[term] != 1 else '') + self.tupleToStr(term)
    return ans
    
  #varsToPlug must be the correct length
  def evaluate(self, varsToPlug):
    s = 0
    nv = len(self.vars)
    for x in self.data:
      term = 1
      for i in xrange(nv):
        term *= varsToPlug[i]**(x[i])
      s += self.data[x] * term
    return s
    
  """  
  def mp_evaluate(self, varsToPlug):  
    s = mpmath.mpf(0)
    nv = len(self.vars)
    for x in self.data:
      term = mpmath.mpf(1)
      for i in xrange(nv):
        term *= mpmath.mpf(varsToPlug[i])**(x[i])
      s += self.data[x] * term
    return s 
  """
       
  def evaluate_symbols(self, varsToPlug, matlab=False):
    s = str(self)
    #s = s.replace(')(', ') (')
    for i, v in enumerate(self.vars):
      s = s.replace('('+v+')', varsToPlug[i])
    return s
  
  def var_degree(self, var):
    return max([t[-1] for t in self.data.keys() if self.data[t] != 0])
    
  def evaluate_polys(self, polys):
    ans = poly('x') ** 0
    ans.data[(0,)] = 0  # make the polynomial "0"
    for tup in self.data:
      tupPoly = poly('x') **0
      tupPoly.data[(0,)] = self.data[tup] # make the multiplier
      for i in xrange(len(tup)):
        tupPoly = tupPoly * (polys[i] ** tup[i])
      ans = ans + tupPoly
    return ans
  
  def deg(self, var=None):
    if var==None:
      return max([sum(k) for k in self.data])
    else:
      ind = self.vars.index(var)
      return max([k[ind] for k in self.data])
    
def distinct_letters(w):
  lw = len(w)
  if lw == 0 or lw == 1:
    return True
  copyw = w.lower()
  for i in xrange(lw-1):
    if copyw[i] in copyw[i+1:]:
      return False
  return True
    
def rotate_left(w, n):
  m = n % len(w)
  return w[m:] + w[:m]
  
def rotate_to_lowest(w):
  ords = map(ord, list(w))
  i = ords.index(min(ords))
  return rotate_left(w, i)
  
def find_repeated_letter_pair(w):
  lw = len(w)
  for i in xrange(lw):
    for j in xrange(i+1, lw):
      if w[i] == w[j]:
        return [i,j]
  return None
  
def find_inverse_letter_pair(w):
  lw = len(w)
  for i in xrange(lw):
    for j in xrange(i, lw):
      if w[i] == w[j].swapcase():
        return [i,j]
  return None 
  

    
#t(x)t(y) = t(xy) + t(xY)  
#returns a polynomial is     
def word_to_trace_poly(w, data={}):
  #print 'Called with ' + w
  if w in data:
    return data[w]
    
  if w == '':
    a = poly('x')
    a.degree = 0
    a.data[tuple([0])] = 2
    return a    
    
  #all distinct letter, but maybe some uppercase
  if distinct_letters(w):
    #if there's an uppercase letter, we want to get rid of it
    if w.islower():
      rw = rotate_to_lowest(w)
      a = poly(rw)
      a.degree = 1
      a.data[tuple([1])] = 1
      return a
    if len(w) == 1:
      rw = w.lower() # tr(a) = tr(A)
      a = poly(rw)
      a.degree = 1
      a.data[tuple([1])] = 1
      return a
    #there must be an uppercase - rotate it to the end, and get rid of it
    i = 0
    while w[i].islower():
      i+=1
    rw = rotate_left(w, i+1)
    w1 = rw[:-1]
    w2 = rw[-1]
    w1w2i = multiply_words(w1, inverse(w2))
  
    w1p = word_to_trace_poly(w1, data)
    w2p = word_to_trace_poly(w2, data)
    w1w2ip = word_to_trace_poly(w1w2i, data)
    
    ans = (w1p*w2p) - w1w2ip
    data[w] = copy.deepcopy(ans)
    return ans
    
  
  #if the letters aren't distinct, find if there are any repeated letters (not inverses)
  t = find_repeated_letter_pair(w)
  if t != None:
    rw = rotate_left(w, t[0]+1)
    t[1] = t[1] - t[0] - 1
    t[0] = len(w)-1
    w1 = rw[:t[1]+1]
    w2 = rw[t[1]+1:]
    w1w2i = multiply_words(w1, inverse(w2))
    
    #print 'rw = ' + rw + '  t = ' + str(t)
    #print 'Calling with ' + w1 + ', ' + w2 + ', ' + w1w2i
    
    w1p = word_to_trace_poly(w1, data)
    w2p = word_to_trace_poly(w2, data)
    w1w2ip = word_to_trace_poly(w1w2i, data)
    
    ans = (w1p*w2p) - w1w2ip 
    
    data[w] = copy.deepcopy(ans)
    
    return ans
    
  #thus there must be a letter and its inverse
  #just make one of them an inverse (choose the uppercase to kill)
  t = find_inverse_letter_pair(w)
  if w[t[0]].isupper():
    rw = rotate_left(w, t[0]+1)
  else:
    rw = rotate_left(w, t[1]+1)
  
  w1 = rw[:-1]
  w2 = rw[-1]
  w1w2i = multiply_words(w1, inverse(w2))
  
  w1p = word_to_trace_poly(w1, data)
  w2p = word_to_trace_poly(w2, data)
  w1w2ip = word_to_trace_poly(w1w2i, data)
  
  ans = (w1p*w2p) - w1w2ip
  
  data[w] = copy.deepcopy(ans)

  return ans
    
    

    
    
    
