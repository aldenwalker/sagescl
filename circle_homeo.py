import copy
import cyclic_order
import string

from sage.all import flatten

class RealInterval:
  """a real interval is just (a,b)"""
  def __init__(self, a, b):
    self.a = a
    self.b = b
    self.size = b-a
  
  def __repr__(self):
    return 'RealInterval(' + str(self.a) + ',' + str(self.b) + ')'
  
  def __str__(self):
    return '[' + str(self.a) + ',' + str(self.b) + ']'
  
  def __contains__(self, x):
    return (self.a <= x) and (x <= self.b)
  
  def get_fraction(self, x):
    return (x-self.a)/self.size
  
  def at_fraction(self, x):
    return self.a + x*self.size


class IntervalMap:
  """This class gives a piecewise linear map between intervals.  It's up 
  to the user to make sure the function is called in its domain/range"""
  def __init__(self, sources=None, targets=None):
    self.sources = copy.deepcopy(sources)
    self.targets = copy.deepcopy(targets)
    self.nints = (len(sources) if sources != None else 0)
  
  def __repr__(self):
    return 'IntervalMap(' + str(self.sources) + ',' + str(self.targets) + ')'
  
  def __str__(self):
    return 'Interval map: ' + ', '.join([str(self.sources[i]) + '->' + str(self.targets[i]) for i in xrange(self.nints)])
  
  def ap(self, x):
    for i in xrange(self.nints):
      if x in self.sources[i]:
        f = self.sources[i].get_fraction(x)
        return self.targets[i].at_fraction(f)
  
  def inverse_ap(self, x):
    for i in xrange(self.nints):
      if x in self.targets[i]:
        f = self.targets[i].get_fraction(x)
        return self.sources[i].at_fraction(f)
      
  def inverse(self):
    return IntervalMap(self.targets, self.sources)

class LinearHomotopyIntervalMap(IntervalMap):
  def __init__(self, I1, I2, t):
    self.I1 = I1
    self.I2 = I2
    self.t = t
    
  def __repr__(self):
    return 'LinearHomotopyIntervalMap('

class EquivariantRHomeo:
  """A real homeo f equivariant under t |-> t+1 is represented 
  by an offset o and a map m:[0,1]->[0,1], where o = f(0).  Then 
  m(x) = f(x) - o, so we have f(x) = f(x/1 + x%1) = m(x%1) + o + x/1.  
  Furthermore, f^-1(y) = (y-o)/1 + m^-1((y-o)%1) """
  
  def __init__(self, offset=None, sources=None, targets=None, imap=None, points=None, nints=None):
    if offset != None:
      self.offset = offset
      if sources != None:
        self.imap = IntervalMap(sources, targets)
      else:
        self.imap = imap
    elif points != None:
      LP = len(points)
      step = 1.0 / LP
      increasing_points = points + [points[0]]
      for i in xrange(LP):
        while increasing_points[i] > increasing_points[i+1]:
          increasing_points[i+1] += 1.0
        while increasing_points[i] < increasing_points[i+1] - 1.0:
          increasing_points[i+1] -= 1.0
      self.offset = increasing_points[0]
      sources = [RealInterval(i*step, (i+1)*step) for i in xrange(LP)]
      targets = [RealInterval(increasing_points[i] - self.offset,   \
                              increasing_points[i+1] - self.offset) for i in xrange(LP)]
      self.imap = IntervalMap(sources, targets)
      
    elif nints != None:
      step = 1.0 / nints
      sources = [RealInterval(i*step, (i+1)*step) for i in xrange(nints)]
      targets = [RealInterval(i*step, (i+1)*step) for i in xrange(nints)]
      self.imap = IntervalMap(sources, targets)
      self.offset = 0
    else:
      self.offset = None
      self.imap = IntervalMap()
  
  def __repr__(self):
    return 'EquivariantRHomeo(offset=' + str(self.offset) + \
                              ', sources=' + str(self.imap.sources) + \
                              ', targets=' + str(self.imap.targets) + ')'
  
  def __str__(self):
    return 'EquivariantRHomeo with offset ' + str(self.offset) + ' and map ' + str(self.imap)
  
  def ap(self, x):
    q, r = divmod(float(x), 1)
    mr = self.imap.ap(r)
    return float(mr + self.offset + q)
  
  def __call__(self, x):
    return self.ap(x)
  
  def inverse_ap(self, x):
    q, r = divmod(float(x)-self.offset, 1)
    mr = self.imap.inverse_ap(r)
    return q + mr
  
  def circle_ap(self, x):
    return self.ap(x) % 1
  
  def inverse_circle_ap(self, x):
    return self.inverse_ap(x) % 1
  
  def inverse(self):
    return SignedEquivariantRHomeo(self, -1)
  
  def __mul__(self, other):
    return ProductEquivariantRHomeo([self, other])
  
  
class SignedEquivariantRHomeo(EquivariantRHomeo):
  def __init__(self, H, s):
    self.H = H
    self.s = s
  
  def __repr__(self):
    return 'SignedEquivariantRHomeo(' + repr(self.H) + ',' + str(self.s) + ')'
  
  def __str__(self):
    return ('+' if self.s>0 else '-') + str(self.H)
  
  def ap(self, x):
    return (self.H.ap(x) if self.s > 0 else self.H.inverse_ap(x))
  
  def __call__(self, x):
    return self.ap(x)
  
  def inverse_ap(self, x):
    return (self.H.inverse_ap(x) if self.s > 0 else self.H.ap(x))
  
  def circle_ap(self, x):
    return (self.H.circle_ap(x) if self.s > 0 else self.H.inverse_circle_ap(x))
  
  def inverse_circle_ap(self, x):
    return (self.H.inverse_circle_ap(x) if self.s > 0 else self.H.circle_ap(x))
  
  def inverse(self):
    return SignedEquivariantRHomeo(self.H, -self.s)
  
class ProductEquivariantRHomeo(EquivariantRHomeo):
  def __init__(self, L):
    if isinstance(L, list):
      all_homeos = [(H.L if isinstance(H, ProductEquivariantRHomeo) else H) for H in L]
      self.L = flatten(all_homeos)
    else:
      self.L = [L]
  
  def __repr__(self):
    return 'Product of ' + str(len(self.L)) + ' homeos'
  
  def __str__(self):
    s = 'Product of '+ str(len(self.L)) + ' homeos:\n'
    for H in self.L:
      s += repr(H) + '\n'
    return s
  
  def ap(self, x):
    ans = x
    for H in reversed(self.L):
      ans = H(ans)
    return ans
  
  def __call__(self, x):
    return self.ap(x)
  
  def inverse_ap(self, x):
    ans = x
    for H in self.L:
      ans = H.inverse_ap(ans)
    return ans
  
  def circle_ap(self, x):
    ans = x
    for H in reversed(self.L):
      ans = H.circle_ap(ans)
    return ans
  
  def inverse_circle_ap(self, x):
    ans = x
    for H in self.L:
      ans = H.inverse_circle_ap(ans)
    return ans

class LinearHomotopyEquivariantRHomeo(EquivariantRHomeo):
  def __init__(self, H0, H1, t):
    self.H0 = H0
    self.H1 = H1
    self.t = t
  
  def __repr__(self):
    return 'LinearHomotopyEquivariantRHomeo(' + repr(self.H0) + ',' + repr(self.H1) + ',' + str(self.t) + ')'
  
  def ap(self, x):
    return (1-t)*self.H0.ap(x) + t*self.H1.ap(x)
  
  def __call__(self, x):
    return self.ap(x)


def PSL2R_action(CO):
  """Gives a list of homeos corresponding to a PSL2R action on the circle
  for the cyclic order given.  All homeos are hyperbolic."""
  
  num_intervals = 2*len(CO)
  gens = string.ascii_lowercase[:len(CO)/2]
  
  #each action takes the source interval to the 
  #complement of the target interval, and puts everything else inside
  #the target interval
  
  homeo_point_maps = []
  step = 1.0/num_intervals
  
  for i in xrange(len(gens)):
    source_i = CO.index(gens[i].swapcase())
    target_i = CO.index(gens[i])
    target_start = target_i*step - step/2.0
    target_end = target_i*step + step/2.0
    target_step = step/(num_intervals-1)
    point_map = [-1 for _ in xrange(num_intervals)]
    for j in xrange(1,num_intervals):
      point_map[(source_i + j)%num_intervals] = target_start + (j-1)*target_step + target_step/2.0
      point_map[source_i] = source_i*step
    homeo_point_maps.append(point_map)
    
  return [EquivariantRHomeo(points=pm) for pm in homeo_point_maps]

      
















  