
import word


class CyclicOrder:
  def __init__(self, inp):
    if type(inp) == str:
      self.CO = inp
    else:
      self.CO = ''.join(inp)
    self.L = len(self.CO)
    self.indices = dict([(g,i) for (i,g) in enumerate(self.CO)])
    
  
  def __repr__(self):
    return "CyclicOrder(" + self.CO + ")"
  
  def __str__(self):
    return self.CO
  
  def __getitem__(self, n):
    if type(n) == slice:
      if type(n.start) == str:
        start = self.indices[n.start]
        stop = self.indices[n.stop]
      else:
        start = n.start
        stop = n.stop
      step = n.step
      if stop < start:
        return self.CO[start::step] + self.CO[:stop:step]
      else:
        return self.CO[start:stop:step]
    return self.CO[n%self.L]
  
  def __contains__(self, x):
    return x in self.indices
  
  def insert(self, location, subslice):
    if type(location) == str:
      loc = self.indices[location]
    else:
      loc = location
    self.__init__(self.CO[:loc+1] + subslice + self.CO[loc+1:])
   
  def __call__(self, a, b, c):
    ai = self.indices.get(a, None)
    if ai == None:
      return 0
    bi = self.indices.get(b, None)
    if bi == None:
      return 0
    ci = self.indices.get(c, None)
    if ci == None:
      return 0
    if bi < ai:
      bi += self.L
    if ci < ai:
      ci += self.L
    return (1 if bi < ci else -1)
  
  def __len__(self):
    return self.L
  
  def step(self, a,b,c):
    """return the number of steps from a to b, without crossing c"""
    ai = self.indices[a]
    bi = self.indices[b]
    ci = self.indices[c]
    if ai < ci:
      ai += self.L
    if bi < ci:
      bi += self.L
    return bi - ai
  
  def rot(self, w):
    if type(w) == list:
      return sum([self.rot(x) for x in w])
    wr = w[::-1] + w[-1]
    s = 0
    for i in xrange(len(w)):
      s += self.step(wr[i], wr[i+1], wr[i+1].swapcase())
    if s%self.L != 0:
      print "Rot isn't an integer?"
      return 0
    return s/self.L
    
  
  
  def extend_to_full_gen_order(self, rank):
    """extend self to a full order on all generators"""
    gens = word.alphabet[:rank]
    gens += [x.swapcase() for x in gens]
    gens_to_add = [g for g in gens if g not in self]
    gens_to_add = ''.join(gens_to_add)
    self.__init__(self.CO + gens_to_add)



def multiple_cyclic_order_eval(t, CO_list) :
  t0, t1, t2 = t
  ords = set([O(t0, t1, t2) for O in CO_list])
  #print "Got orders ", ords, " for tripod ", (t0, t1, t2), " under orders: ", CO_list
  if 1 in ords:
    if -1 in ords:
      return None
    else: 
      return 1
  elif -1 in ords:
    return -1
  return 0

def four_tuple_from_cyclic_orders(t, CO_list):
  s = [ multiple_cyclic_order_eval(t[:j] + t[j+1:], CO_list) for j in xrange(0, 4)]
  #print "I evaluated the 4-tuple ", t, " on the cyclic orders and got: ", str(s)
  if None in s:
    return None
  got_order = None
  #0123 and reverse
  if s[0] == s[2] and s[0] != 0:
    got_order = (CyclicOrder(t) if s[0] == 1 else CyclicOrder(t[::-1]))
  elif s[1] == s[3] and s[1] != 0:
    got_order = (CyclicOrder(t) if s[1] == 1 else CyclicOrder(t[::-1]))
  #0132 and reverse
  elif s[2] == -s[1] and s[2] != 0:
    got_order = CyclicOrder( [t[j] for j in ([0,1,3,2] if s[2]==1 else [2,3,1,0])] )
  elif s[3] == -s[0] and s[3] != 0:
    got_order = CyclicOrder( [t[j] for j in ([0,1,3,2] if s[3]==1 else [2,3,1,0])] )
  #0312 and reverse
  elif s[0] == -s[1] and s[0] != 0:
    got_order = CyclicOrder( [t[j] for j in ([0,3,1,2] if s[0]==1 else [2,1,3,0])] )
  elif s[3] == -s[2] and s[3] != 0:
    got_order = CyclicOrder( [t[j] for j in ([0,3,1,2] if s[3]==1 else [2,1,3,0])] )
  
  #print "I got the order: ", str(got_order)
  
  #we didn't conclude anything
  if got_order == None:
    return False

  #check consistency
  s2 = [ got_order(*(t[:j] + t[j+1:])) for j in xrange(0, 4) ]
  for j in xrange(4):
    if s[j] != 0 and s2[j] != s[j]:
      return None
  return got_order

def sorted_partial_order(L, f):
  """sorts L relative to the function f, where 
  f(x,y) < 0 if x < y, ==0 if equal or unknown, and 
  f(x,y) > 0 if x > y.  it basically has to be insertion sort"""
  LL = list(L)
  for i in xrange(1,len(L)):
    j=i-1
    while j >= 0 and f(LL[j], LL[j+1]) >= 0:
      temp = LL[j+1]
      LL[j+1] = LL[j]
      LL[j] = temp
      j -= 1
  return LL
    

def extend_suborders_to_order(rank, T):
  #print "Called to extend the orders: ", T
  gens = word.alphabet[:rank]
  gens += [x.swapcase() for x in gens]
  def cmp_rel_gen0(g1, g2):
    return multiple_cyclic_order_eval([gens[0], g2, g1], T)
  gens = [gens[0]] + sorted_partial_order(gens[1:], cmp_rel_gen0)
  #print "Sorted order to ", gens
  return CyclicOrder(gens)
