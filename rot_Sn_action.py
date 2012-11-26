from sage.all import *

from word import *


alphabet = list('abcdefghijklmnopqrstuvwxyz')

def all_cyclic_orders(rank):
  letters = alphabet[:rank] + inverse(alphabet[:rank])
  P = permutations(letters[1:])
  return [''.join([letters[0]] + p) for p in P]

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

def tripod_boundary(t, all_boundary_words=None):
  if all_boundary_words==None:
    return [inverse(t[0]) + t[1], inverse(t[1]) + t[2], inverse(t[2]) + t[0]]
  else:
    all_bds = []
    for b in [inverse(t[0]) + t[1], inverse(t[1]) + t[2], inverse(t[2]) + t[0]]:
      for i in xrange(len(b)-all_boundary_words+1):
        all_bds.append( b[i:i+all_boundary_words] )
    return all_bds

def tripods(rank, edge_len=1):
  gens = alphabet[:rank] + inverse(alphabet[:rank])
  unordered = sage.combinat.subset.Subsets_sk(gens, 3)
  tripods_len_1 = [list(t) for t in unordered] + [[t[0], t[2], t[1]] for t in unordered]
  if edge_len==1:
    return tripods_len_1
  #otherwise, the set of all tripods is all possible
  #extensions of these tripods of length 1
  W = all_words_of_len(edge_len, alphabet[:rank])
  starting_with = dict( [ (g, [w for w in W if w[0]==g]) for g in alphabet[:rank] + inverse(alphabet[:rank])] )
  tripods_len_ell = []
  for T in tripods_len_1:
    word_lists = [ starting_with[T[0]], starting_with[T[1]], starting_with[T[2]] ]
    word_list_lens = len(word_lists[0]) #they should all have the same length
    selected = tuples(range(word_list_lens), 3)
    tripods_len_ell += [ [word_lists[0][s[0]], word_lists[1][s[1]], word_lists[2][s[2]]] for s in selected]
  return tripods_len_ell
    

#word_ind_dict should be a dictionary that maps each word to 
#+/i, where +i if the word lives in index i-1 and -i if the word 
#lives in i-1 but with negative sign
def words_to_vector(C, L, word_ind_dict):
  ans = [0 for i in xrange(L)]
  for w in C:
    widw = word_ind_dict[w]
    if widw < 0:
      ans[-widw-1] -= 1
    else:
      ans[widw-1] += 1
  return ans

class CQ:
  @staticmethod
  def unit_ball_len_2(rank, W_in=None, back='ppl'):
    if W_in == None:
      W = CQ.antisym_basis(2, rank)
    else:
      W = W_in
    if (2,rank) not in CQ.homs:
      CQ.compute_rank_homs(2, rank)
    equalities_list = [[0] + c.to_vector(W) for c in CQ.homs[(2,rank)]]
    T = tripods(rank)
    word_ind_dict = {}
    for i in xrange(len(W)):
      word_ind_dict[W[i]] = i+1
      word_ind_dict[inverse(W[i])] = -(i+1)
    ineqs_list = []
    ambient_dim = len(W)
    for t in T:
      bd = tripod_boundary(t)
      vec = words_to_vector(bd, ambient_dim, word_ind_dict)
      negvec = [-v for v in vec]
      ineqs_list.append([1] + negvec)
    r = (QQ if back=='ppl' or back=='cddr' else RDF)
    UB = Polyhedron(base_ring=r, ieqs=ineqs_list, eqns=equalities_list, backend=back)
    return UB

  @staticmethod
  def unit_ball(ell, rank, W_in=None, back='ppl'):
    if W_in == None:
      W = CQ.antisym_basis(ell, rank)
    else:
      W = W_in
    if (ell,rank) not in CQ.homs:
      CQ.compute_rank_homs(ell, rank)
    hom_equalities_list = [[0] + c.to_vector(W) for c in CQ.homs[(ell,rank)]]
    if (ell, rank) not in CQ.zero_vectors:
      CQ.compute_zero_vectors(ell, rank)
    zero_equalities_list = [[0] + c.to_vector(W) for c in CQ.zero_vectors[(ell,rank)]]
    T = tripods(rank, ell-1)
    word_ind_dict = {}
    for i in xrange(len(W)):
      word_ind_dict[W[i]] = i+1
      word_ind_dict[inverse(W[i])] = -(i+1)
    ineqs_list = []
    ambient_dim = len(W)
    for t in T:
      bd = tripod_boundary(t, all_boundary_words=ell)
      vec = words_to_vector(bd, ambient_dim, word_ind_dict)
      negvec = [-v for v in vec]
      ineqs_list.append([1] + negvec)
    r = (QQ if back=='ppl' or back=='cddr' else RDF)
    UB = Polyhedron(base_ring=r, ieqs=ineqs_list, eqns=hom_equalities_list + zero_equalities_list, backend=back)
    return UB

  @staticmethod
  def from_vector(ell, rank, v, W):
    d = dict( [ (W[i], v[i]) for i in xrange(len(W))] )
    return CQ(ell=ell, rank=rank, rules=d)
    

  @staticmethod
  def antisym_basis(ell, rank):
    W = all_words_of_len(ell, alphabet[:rank])
    pairs = [(x, inverse(x)) for x in W]
    pairs = list(set([ tuple(sorted(p, reverse=True)) for p in pairs]))
    return [p[0] for p in pairs]

  @staticmethod
  def compute_rank_homs(ell, rank):
    if (ell, rank) in CQ.homs:
      return
    CQ.homs[(ell, rank)] = []
    W = []
    for i in xrange(ell):
      W.append(all_words_of_len(i, list(alphabet[:rank])))
    for g in alphabet[:rank]:
      for i in xrange(ell):
        G = g.swapcase()
        current_words = []
        for w1 in W[i]:
          if len(w1) > 0 and w1[-1] == G:
            continue
          for w2 in W[ell-i-1]:
            if len(w2) > 0 and w2[0] == G:
              continue
            current_words.append(w1+g+w2)
        #print "Current words for ", g, " for i=", i,": ", current_words
        d = dict( [(w,1) for w in current_words] )
        CQ.homs[(ell,rank)].append(CQ(ell=ell, rank=rank, rules=d))

  @staticmethod
  def compute_zero_vectors(ell, rank):
    #produces CQs for each word of length ell-1, which is (everything ending with it) - 
    #(everything beginning with it).  These are zero on all chains
    CQ.zero_vectors[(ell, rank)] = []
    Wmo = all_words_of_len(ell-1, alphabet[:rank])
    gens = alphabet[:rank] + inverse(alphabet[:rank])
    for w in Wmo:
      all_beginning = [w + g for g in gens if w[-1] != g.swapcase()]
      all_ending = [g + w for g in gens if w[0] != g.swapcase()]
      pairs = [(e,1) for e in all_ending] + [(b,-1) for b in all_beginning]
      CQ.zero_vectors[(ell,rank)].append( CQ(ell=ell, rank=rank, rules=dict(pairs)) )


  homs = {}
  zero_vectors = {}
  
  def __init__(self, order=None, ell=None, rank=None, rules=None):
    if order == None:
      if rules == None:
        self.rules = {}
        self.ell = ell
        self.rank = rank
      else:
        self.rules = rules
        self.ell = ell
        self.rank = rank
        self.antisymmetrize()
      return
    #otherwise, we're creating it from a cyclic order
    self.rules = {}
    self.rank = chain_rank([order])
    self.ell = 2
    allgens = alphabet[:self.rank] + inverse(alphabet[:self.rank])
    for g1 in allgens:
      g1s = g1.swapcase()
      for g2 in allgens:
        if g2 == g1 or g2 == g1s:
          continue
        self.rules[g1+g2] = cyclic_order_coefficient(order, g1, g2, g2.swapcase())
    self.antisymmetrize()


  def __iadd__(self, other):
    for x in other.rules:
      if x in self.rules:
        self.rules[x] += other.rules[x]
      elif inverse(x) in self.rules:
        self.rules[inverse(x)] -= other.rules[x]
    return self

  def __add__(self, other):
    d = {}
    all_keys = list(set(self.rules.keys() + other.rules.keys()))
    for k in all_keys:
      d[k] = self.rules.get(k, 0) + other.rules.get(k, 0)
    return CQ(ell=self.ell, rank=self.rank, rules=d)

  def __rmul__(self, other):
    d = {}
    for x in self.rules:
      d[x] = self.rules[x] * other
    return CQ(ell=self.ell, rank=self.rank, rules=d)
  
  def __repr__(self):
    return "CQ(ell="+str(self.ell)+", rank="+str(self.rank)+", rules="+str(self.rules)+")"

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
  
  def increase_ell(self, new_ell):
    new_rules = []
    W = all_words_of_len(new_ell-self.ell, alphabet[:self.rank])
    for w in self.rules:
      new_rules += [ (w+w2, self.rules[w]) for w2 in W if w2[0] != w[-1].swapcase()] 
    return CQ(ell=new_ell, rank=self.rank, rules=dict(new_rules))
    

  def add_words(self, d):
    for x in d:
      self.rules[x] = self.rules.get(x,0) + d[x]
    self.antisymmetrize()

  def antisymmetrize(self):
    #first antisymmetrize
    todo = self.rules.keys()
    pairs = [(x, inverse(x)) for x in todo]
    pairs = [tuple(sorted(p, reverse=True)) for p in pairs]
    pairs = list(set(pairs))
    for p in pairs:
      coef = self.rules.get(p[0], 0) - self.rules.get(p[1],0)
      if coef == 0:
        if p[0] in self.rules:
          del self.rules[p[0]]
      else:
        self.rules[p[0]] = coef
      if p[1] in self.rules:
        del self.rules[p[1]]

  #we need to determine if their difference is in the span of the homs
  def __eq__(self, other, W_in=None):
    if (self.ell, self.rank) not in CQ.homs:
      CQ.compute_rank_homs(self.ell, self.rank)
    dif = -1*other
    dif = self + dif
    #print "Difference is: ", dif
    if W_in == None:
      W = CQ.antisym_basis(self.ell, self.rank)
    else:
      W = W_in
    hom_basis = [c.to_vector(W) for c in CQ.homs[(self.ell, self.rank)]]
    my_vector = dif.to_vector(W)
    A = Matrix(QQ, hom_basis)
    #print "Matrix: "
    #print str(A)
    #print "diff vector: "
    #print str(my_vector)
    try:
      X = A.solve_left(vector(my_vector))
      return True
    except ValueError:
      return False
    
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
      
  def __str__(self):
    return str(self.rules)
  
  def ap(self, O):
    return ''.join([self.rules[g] for g in O])
  


def find_rots_spanning(rank):
  lets = alphabet[:rank] + inverse(alphabet[:rank])
  M = []
  good_rots = []
  W = CQ.antisym_basis(2, rank)
  if (2,rank) not in CQ.homs:
    CQ.compute_rank_homs(2, rank)
  hom_basis = [c.to_vector(W) for c in CQ.homs[(2,rank)]]
  V = VectorSpace(QQ, len(W))
  HS = V.subspace([ V(h) for h in hom_basis])
  QH = V.quotient(HS)
  target_rank = QH.dimension()
  while True:
    orders = []
    for i in xrange(target_rank):
      lets2 = lets
      shuffle(lets2)
      lets2 = ''.join(lets2)
      orders.append(lets2)
    rots = [CQ(O) for O in orders]
    M = Matrix(QQ, [QH(r.to_vector(W)) for r in rots])
    if M.rank() == target_rank:
      break
  return (V, QH, orders, rots, M)

def is_S2k_linear(rank):
  V, QH, orders, rots, M = find_rots_spanning(rank)
  perm1 = ''.join([alphabet[i] + inverse(alphabet[i]) for i in xrange(rank)])
  perm2 = 'Aa' + perm1[2:]
  perm1 = 'bAaB' + perm1[4:]
  p1 = CO_Perm(perm1)
  p2 = CO_Perm(perm2)
  W = CQ.antisym_basis(2, rank)
  all_O = all_cyclic_orders(rank)
  for p in (p1, p2):
    target_rots = [CQ(p.ap(O)) for O in orders]
    target_matrix = Matrix(QQ, [QH(r.to_vector(W)) for r in target_rots])
    transition_matrix = M.inverse() * target_matrix
    for O in all_O:
      r1 = CQ(O)
      r2 = CQ(p.ap(O))
      v1 = QH(r1.to_vector(W))
      v2 = QH(r2.to_vector(W))
      if v1 * transition_matrix != v2:
        print "I think that perm ", p, " doesn't work with rot ", O
        print 1/0
        return False
  return True
              
def rots_convex_hull(rank, br=QQ, back='ppl'):
  rots = [CQ(O) for O in all_cyclic_orders(rank)]
  W = CQ.antisym_basis(2, rank)
  rots_vecs = [r.to_vector(W) for r in rots]
  if (2,rank) not in CQ.homs:
    CQ.compute_rank_homs(2, rank)
  H = [h.to_vector(W) for h in CQ.homs[(2,rank)]]
  V = VectorSpace(QQ, len(W))
  HS = V.subspace(H)
  QH = V.quotient(HS)
  rots_quotient = [QH(r) for r in rots_vecs]
  return Polyhedron(base_ring=br, backend=back, vertices=rots_quotient)

def save_rots_convex_hull(rank, filename, br, back):
  rots_convex_hull(rank, br, back).save(filename)


















