import copy
import cyclic_order
import string
import math


def dedup(L):
  ans = []
  for x in L:
    if ans == [] or abs(ans[-1]-x) > 0.0000001:
      ans.append(x)
  return ans
  
def merge(L1, L2):
  ans = []
  i1 = 0
  i2 = 0
  LL1 = len(L1)
  LL2 = len(L2)
  while True:
    if i1 >= LL1:
      for i in xrange(i2, LL2):
        if ans == [] or abs(ans[-1]-L2[i]) > 0.0000001:
          ans.append(L2[i])
      break
    elif i2 >= LL2:
      for i in xrange(i1, LL1):
        if ans == [] or abs(ans[-1]-L1[i]) > 0.0000001:
          ans.append(L1[i])
      break
    else:
      if L1[i1] < L2[i2]:
        x = L1[i1]
        i1 += 1
      else:
        x = L2[i2]
        i2 += 1
      if ans == [] or abs(ans[-1]-x) > 0.0000001:
        ans.append(x)
  return ans

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
    eqa = (abs(self.a-x) < 0.0000001)
    eqb = (abs(self.b-x) < 0.0000001)
    return (eqa or self.a <= x) and (eqb or x <= self.b)
  
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
    
  def breakpoints(self):
    return [I.a for I in self.sources] + [self.sources[-1].b]

def IM_homotopy(m1, m2, t):
  """returns the IntervalMap which is (1-t)*m1 + t*m2"""
  all_breaks = merge(m1.breakpoints(), m2.breakpoints())
  new_sources = [RealInterval(all_breaks[i], all_breaks[i+1]) for i in xrange(len(all_breaks)-1)]
  new_images = [(1-t)*m1.ap(x) + t*m2.ap(x) for x in all_breaks]
  new_targets = [RealInterval(new_images[i], new_images[i+1]) for i in xrange(len(all_breaks)-1)]
  return IntervalMap(sources=new_sources, targets = new_targets)

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
    all_breaks = [0,1] + [self.ap(x)%1 for x in self.imap.breakpoints()]
    all_breaks.sort()
    all_breaks = dedup(all_breaks)
    nbreaks = len(all_breaks)
    new_offset = self.inverse_ap(0)
    new_sources = [RealInterval(all_breaks[i], all_breaks[i+1]) for i in xrange(nbreaks-1)]
    new_images = [self.inverse_ap(x)-new_offset for x in all_breaks]
    new_targets = [RealInterval(new_images[i], new_images[i+1]) for i in xrange(nbreaks-1)]
    return EquivariantRHomeo(offset=new_offset, sources=new_sources, targets=new_targets)
    
  
  def __mul__(self, other):
    """returns self*other; note (self*other)(x) = self(other(x))"""
    pullback_breakpoints = [other.inverse_ap(x)%1 for x in self.imap.breakpoints()]
    pullback_breakpoints.sort()
    #print "Merging breakpoints ", other.imap.breakpoints(), pullback_breakpoints
    all_breakpoints = merge(other.imap.breakpoints(), pullback_breakpoints)
    npoints = len(all_breakpoints)
    new_offset = self.ap(other.ap(0))
    new_images = [self.ap(other.ap(x))-new_offset for x in all_breakpoints]
    #print "new breakpoints: ", all_breakpoints
    #print "new images: ", new_images
    new_sources = [RealInterval(all_breakpoints[i], all_breakpoints[i+1]) for i in xrange(npoints-1)]
    new_targets = [RealInterval(new_images[i], new_images[i+1]) for i in xrange(npoints-1)]
    return EquivariantRHomeo(offset=new_offset, sources=new_sources, targets=new_targets)
  
  def rot(self, n=100):
    x = 0
    for i in xrange(n):
      x = self.ap(x)
    return x/100.0

def ERH_homotopy(h1, h2, t):
  """returns the map which is (1-t)*h1 + t*h2"""
  new_imap = IM_homotopy(h1.imap, h2.imap, t)
  new_offset = (1-t)*h1.offset + t*h2.offset
  return EquivariantRHomeo(offset=new_offset, imap=new_imap)

def L2_distance(h1, h2):
  all_breaks = merge(h1.imap.breakpoints(), h2.imap.breakpoints())
  s = 0
  for i in xrange(len(all_breaks)-1):
    x1 = all_breaks[i]
    x2 = all_breaks[i+1]
    y11 = h1.ap(x1)
    y12 = h2.ap(x1)
    y21 = h1.ap(x2)
    y22 = h2.ap(x2)
    a = (y21-y11)/(x2-x1)
    b = y11 - x1*a
    c = (y22-y12)/(x2-x1)
    d = y12 - x1*c
    if abs(a-c) < 0.0000001:
      s += (b-d)*(b-d)*(x2-x1)
      continue
    I1 = ((a-c)*x1 + (b-d))**3/(3*(a-c))
    I2 = ((a-c)*x2 + (b-d))**3/(3*(a-c))
    s += I2-I1
  return math.sqrt(s)


def homotope_and_compare(H, gen_words, target_word,n):
  """homotope gens between the identity and the words in gen_word, 
  then compute target_word in *those* generators, and compare it to 
  H[0]"""
  ALC = string.ascii_lowercase
  gen_words_L = [[(ALC.index(x)+1 if x.islower() else -(ALC.index(x.lower())+1)) for x in w] for w in gen_words]
  target_word_L = [(ALC.index(x)+1 if x.islower() else -(ALC.index(x.lower())+1)) for x in target_word]
  step = 1.0/n
  h = [[EquivariantRHomeo(nints=8) for __ in xrange(2)] for _ in xrange(2)]
  for j in xrange(2):
    for i in xrange(len(gen_words_L[j])):
      if gen_words_L[j][i] > 0:
        h[j][1] = h[j][1]*H[gen_words_L[j][i]-1]
      else:
        h[j][1] = h[j][1]*H[-gen_words_L[j][i]-1].inverse()
  
  current_h = [0,0]
  plot_data = []
  plot_data_rot = []
  for i in xrange(n+1):
    current_h[0] = ERH_homotopy(h[0][0], h[0][1], i*step)
    for j in xrange(n+1):
      current_h[1] = ERH_homotopy(h[1][0], h[1][1], j*step)
      current_prod = EquivariantRHomeo(nints=8)
      for k in xrange(len(target_word)):
        if target_word_L[k] > 0:
          current_prod = current_prod*current_h[target_word_L[k]-1]
        else:
          current_prod = current_prod*current_h[-target_word_L[k]-1].inverse()
      plot_data.append( (i*step, j*step, L2_distance(current_prod, H[0])) )
      plot_data_rot.append( (i*step, j*step, current_prod.rot()) )
  return (plot_data, plot_data_rot)


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

      
















  
