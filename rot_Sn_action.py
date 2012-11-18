from sage.all import *

from word import *


alphabet = list('abcdefghijklmnopqrstuvwxyz')


def cyclic_order_coefficient(O, a, b, c):
  LO = len(O)
  ai = O.index(a)
  bi = O.index(b)
  ci = O.index(c)
  if ai < ci:
    ai += LO
  if bi < ci:
    bi += LO
  return bi-ai


class CQ:
  homs = {}
  
  def __init__(self, order=None, ell=None, rank=None, rules=None, uniq=True):
    if order == None:
      if rules == None:
        self.rules = {}
        self.ell = ell
        self.rank = rank
      else:
        self.rules = rules
        self.ell = ell
        self.rank = rank
        if uniq:
          self.uniqueify()
      return
    #otherwise, we're creating it from a cyclic order
    self.rules = {}
    self.rank = chain_rank([order])
    self.ell = 2
    for g1 in alphabet[:self.rank]:
      g1s = g1.swapcase()
      for g2 in alphabet[:self.rank]:
        if g2 == g1 or g2 == g1s:
          continue
        self.rules[g1+g2] = cyclic_order_coefficient(order, g1, g2, g2.swapcase())
    self.uniqueify()


  def __iadd__(self, other):
    for x in other.rules:
      if x in self.rules:
        self.rules[x] += other.rules[x]
      elif inverse(x) in self.rules:
        self.rules[inverse(x)] -= other.rules[x]
    return self

  def __rmul__(self, other):
    for x in self.rules:
      self.rules[x] *= other
  
  def __repr__(self):
    return "CQ(ell="+str(self.ell)+"rank="+str(self.rank)+"rules="+str(self.rules)+")"

  def __str__(self):
    return self.__repr__()

  def ap(self, w, kind='chain'):
    if type(w) == list:
      return sum([self.ap(x, kind) for x in w])
    Lw = len(w)
    s = 0
    if kind == 'chain':
      for i in xrange(Lw):
        x = cyclic_subword(w, i, self.ell)
        s += self.rules.get(x, 0)
        s -= self.rules.get(inverse(x), 0)
      return s
    elif kind == 'word':
      for i in xrange(Lw):
        if i+self.ell > Lw:
            break
        x = cyclic_subword(w, i, self.ell)
        s += self.rules.get(x, 0)
        s -= self.rules.get(inverse(x), 0)
      return s
    else:
      print "kind not recognized"
  
  def add_words(self, d):
    for x in d:
      self.rules[x] = self.rules.get(x,0) + d[x]
    self.uniqueify()

  def uniqueify(self):
    #first antisymmetrize
    todo = self.rules.keys()
    pairs = [(x, inverse(x)) for x in todo]
    pairs = [sorted(p, reverse=True) for p in pairs]
    pairs = list(set(pairs))
    for p in todo:
      self.rules[p[0]] = self.rules.get(p[0], 0) - self.rules.get(p[1],0)
    #now remove coefficients of aa, bb, etc
    if self.rank not in CQ.homs:
      self.compute_rank_homs()
    for x in alphabet[:rank]:
      xx = x+x
      if xx not in self.rules:
        continue
      coef = self.rules[xx]
      self += ((-coef) * CQ.homs[self.rank][xx])

  def compute_rank_homs(self):
    if self.rank in CQ.homs:
      return
    CQ.homs[self.rank] = {}
    W = all_words_of_len(self.len-1, list(alphabet[:self.rank]))
    for g in alphabet[:self.rank]:
      G = g.swapcase()
      my_words = [g + w for w in W if w[-1] != G]
      d = dict( [(w, 1) for w in my_words] )
      CQ.homs[self.rank][g+g] = CQ(self.len, self.rank, rules=d, uniq=False)
  
  def to_vector(self, W):
    return [self.rules.get(w,0) for w in W]




#a permutation class to act on cyclic orders
#giving a cyclic order specs the perm by sending aAbBcC... to the order
class CO_Perm:
  def __init__(self, CO=None, rank=None):
    if CO==None:
      if rank==None:
        self.rules=None
        self.rank=None
        return
      self.rank = rank
      self.rules = dict( [(g,g) for g in alphabet[:self.rank]] + \
                         [(g.swapcase(),g.swapcase()) for g in alphabet[:self.rank]])
      return
    else:
    self.rank = chain_rank([CO])
    self.rules = {}
    for i in xrange(self.rank):
      self.rules[alphabet[i]] = CO[2*i]
      self.rules[inverse(alphabet[i])] = CO[2*i+1]
    return
      
  def ap(self, O):
    return ''.join([self.rules[g] for g in O])
  




















