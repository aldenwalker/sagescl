from sage.all import *

from word import *
from morph import *

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

def all_distinct(L):
  return len(set(L)) == len(L)

def tripod_forms_sphere(t):
  return all_distinct([t[0][-1], t[1][-1], t[2][-1]])

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
    
#returns a list of tripods which should all be positive if the quadpod is in order
def tripods_from_quadpod(q):
  return [ [q[0], q[1], q[2]], [q[2],q[3],q[0]], [q[3],q[0],q[1]], [q[1],q[2],q[3]] ]

#given a dictionary which assigns +/- 1 to tripods (0 on the others), 
#say whether it is compatible on quadpods
def compatible_on_quadpods(T_dict_in, rank, Q_in=None, T_in=None):
  gens = alphabet[:rank] + inverse(alphabet[:rank])
  if Q_in == None:
    Q = sage.combinat.subset.Subsets_sk(gens, 4)
  else:
    Q = Q_in
  if T_in == None:
    T = tripods(rank)
  else:
    T = T_in
  T_dict = {}
  for t in T_dict_in:
    T_dict[t] = T_dict_in[t]
    T_dict[(t[0], t[2], t[1])] = -T_dict_in[t]
    T_dict[(t[1], t[2], t[0])] = T_dict_in[t]
    T_dict[(t[1], t[0], t[2])] = -T_dict_in[t]
    T_dict[(t[2], t[0], t[1])] = T_dict_in[t]
    T_dict[(t[2], t[1], t[0])] = -T_dict_in[t]
  for q in Q:
    #for each quadruple, check that it is compatible
    #I'm just going to do this brute force
    all_quadpods = [ [q[0]] + x for x in permutations(q[1:])]
    its_compatible = False
    for qp in all_quadpods:
      tfq = tripods_from_quadpod(qp)
      if all([v==0 or v==1 for v in [T_dict.get(tuple(t), 0) for t in tfq]]):
        its_compatible = True
        break
    if not its_compatible:
      break
  return its_compatible

#returns the four boundary words of length 2
def quadpod_boundary(q):
  return [q[0].swapcase() + q[1], \
          q[1].swapcase() + q[2], \
          q[2].swapcase() + q[3], \
          q[3].swapcase() + q[0]]

#given a dictionary which assigns +/-1 to tripods, 
#this gives the counting function which is the boundary of all the 
#quadpods which are compatible
def quadpod_boundary_function(T_dict_in, rank, Q_in=None, T_in=None):
  gens = alphabet[:rank] + inverse(alphabet[:rank])
  if Q_in == None:
    Q = sage.combinat.subset.Subsets_sk(gens, 4)
  else:
    Q = Q_in
  if T_in == None:
    T = tripods(rank)
  else:
    T = T_in
  T_dict = {}
  for t in T_dict_in:
    T_dict[t] = T_dict_in[t]
    T_dict[(t[0], t[2], t[1])] = -T_dict_in[t]
    T_dict[(t[1], t[2], t[0])] = T_dict_in[t]
    T_dict[(t[1], t[0], t[2])] = -T_dict_in[t]
    T_dict[(t[2], t[0], t[1])] = T_dict_in[t]
    T_dict[(t[2], t[1], t[0])] = -T_dict_in[t]
  all_words = {}
  for q in Q:
    #for each quadruple, check that it is compatible
    #I'm just going to do this brute force
    all_quadpods = [ [q[0]] + x for x in permutations(q[1:])]
    its_compatible = False
    compatible_quadpods = []
    for qp in all_quadpods:
      tfq = tripods_from_quadpod(qp)
      if all([v==0 or v==1 for v in [T_dict.get(tuple(t), 0) for t in tfq]]):
        its_compatible = True
        compatible_quadpods.append(qp)
    if not its_compatible:
      break
    #print "For ", q, " found the compatible quadpods:"
    #print compatible_quadpods
    #add a boundary word wherever it is needed (part of a positive tripod)
    num_cqp = len(compatible_quadpods)
    for cqp in compatible_quadpods:
      #print "For ", cqp, " adding ",
      for i in xrange(4):
        a,b,c,d = cqp[i:] + cqp[:i]
        if True: #T_dict.get((a,b,c), 0) == 1 or T_dict.get((a,b,d),0)==1:
          w = a.swapcase() + b
          all_words[w] = all_words.get(w,0)+(1/Integer(num_cqp))
          #print w, (1/Integer(num_cqp)), 
        else:
          pass
          #print "Not ", a.swapcase() + b, (1/Integer(num_cqp)),

  if not its_compatible:
    return None
  return CQ(ell=2, rank=rank, rules=all_words)
  

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
      neg2vec = [-2*v for v in vec]
      ineqs_list.append([1] + neg2vec)
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
      neg2vec = [-2*v for v in vec]
      ineqs_list.append([1] + neg2vec)
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
  def compute_QH_quotient(ell, rank):
    W = CQ.antisym_basis(ell, rank)
    if (ell, rank) not in CQ.homs:
      CQ.compute_rank_homs(ell, rank)
    if (ell, rank) not in CQ.zero_vectors:
      CQ.compute_zero_vectors(ell, rank)
    Z = [c.to_vector(W) for c in CQ.zero_vectors[(ell, rank)]]
    H = [c.to_vector(W) for c in CQ.homs[(ell, rank)]]
    V = VectorSpace(QQ, len(W))
    ZH = V.subspace(Z+H)
    return W, V.quotient(ZH)

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
      L = [x for x in all_beginning if x in all_ending]
      all_beginning = [x for x in all_beginning if x not in L]
      all_ending = [x for x in all_ending if x not in L]
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

  def defect(self, T_in=None):
    if T_in == None:
      T = tripods(self.rank, self.ell-1)
    else:
      T = T_in
    return 2*max([self.ap(tripod_boundary(t, all_boundary_words=self.ell), kind='word') for t in T])

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
    
  def tripod_ap(self, t):
    return self.ap(tripod_boundary(t, all_boundary_words=self.ell), kind='word')
  
  def compatible_on_quadpods(self, Q_in=None, T_in=None):
    if self.ell > 2:
      print "Not implemented for ell > 2"
      return
    gens = alphabet[:self.rank] + inverse(alphabet[:self.rank])
    if Q_in == None:
      Q = sage.combinat.subset.Subsets_sk(gens, 4)
    else:
      Q = Q_in
    if T_in == None:
      T = tripods(self.rank)
    else:
      T = T_in
    D = self.defect(T)
    for q in Q:
      #for each quadruple, check that it is compatible
      #I'm just going to do this brute force
      all_quadpods = [ [q[0]] + x for x in permutations(q[1:])]
      its_compatible = False
      for qp in all_quadpods:
        tfq = tripods_from_quadpod(qp)
        if all([v==0 or 2*v==D for v in [self.tripod_ap(t) for t in tfq]]):
          its_compatible = True
          break
      if not its_compatible:
        break
    return its_compatible
      

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
    if (self.ell, self.rank) not in CQ.zero_vectors:
      CQ.compute_zero_vectors(self.ell, self.rank)
    dif = -1*other
    dif = self + dif
    #print "Difference is: ", dif
    if W_in == None:
      W = CQ.antisym_basis(self.ell, self.rank)
    else:
      W = W_in
    hom_basis = [c.to_vector(W) for c in CQ.homs[(self.ell, self.rank)]]
    zero_basis = [c.to_vector(W) for c in CQ.zero_vectors[(self.ell, self.rank)]]
    my_vector = dif.to_vector(W)
    A = Matrix(QQ, hom_basis + zero_basis)
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

  def act_by_aut(self, A):
    #we need to figure out how big the tripods need to be 
    #to determine tripods of length ell-1
    #we need at least ell-1, obviously
    current_edge_len_attempt = self.ell-1
    gens = alphabet[:self.rank] + inverse(alphabet[:self.rank])
    while True:
      W = all_words_of_len(current_edge_len_attempt, gens)
      by_first_letter = dict([ (g,[w for w in W if w[0] == g]) for g in gens])
      sig_by_first_letter = dict([ (g, [sig_word_under_aut(w, A, rank=self.rank) for w in by_first_letter[g]]) for g in gens])
      min_leftover = 100*self.ell
      for g1 in gens:
        for w1 in sig_by_first_letter[g1]:
          Lw1 = len(w1)
          for g2 in gens:
            if g2 == g1:
              continue
            for w2 in sig_by_first_letter[g2]:
              Lw2 = len(w2)
              leftover = min(Lw1, Lw2)
              m = leftover
              i=0
              while i < m and w1[i] == w2[i]:
                leftover -= 1
                i += 1
              if leftover < min_leftover:
                min_leftover = leftover
      if min_leftover >= self.ell-1:
        break
      else:
        current_edge_len_attempt += 1
        if current_edge_len_attempt > 2:
          return None
    #if we get here, then current_edge_len_attempt gives us what the 
    #tripods are of length self.ell-1  now we need to actually act by them 
    edge_len = current_edge_len_attempt
    T = tripods(self.rank, edge_len)
    W = CQ.antisym_basis(edge_len+1, self.rank)
    dim = len(W)
    rows = []
    values = []
    for t in T:
      image_t = A.act_on_tripod(t)
      if min(map(len, image_t)) < self.ell-1:
        print "Error -- tripod too small"
      image_t = [w[:self.ell-1] for w in image_t]
      val = self.ap(tripod_boundary(image_t), kind='word')
      word_vec = tripod_boundary(t, all_boundary_words=edge_len+1)
      vec = [0 for i in xrange(dim)]
      for i in xrange(dim):
        vec[i] += word_vec.count(W[i]) - word_vec.count(inverse(W[i]))
      values.append(val)
      rows.append(vec)
    M = Matrix(QQ, rows)
    b = vector(QQ, values)
    try:
      soln = M.solve_right(b)
    except:
      return None
    rules = dict([(W[i], soln[i]) for i in xrange(dim)])
    return CQ(ell=edge_len+1, rank=self.rank, rules=rules)
      
     


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
  
  def ap(self, O, n=1):
    ans = O
    for i in xrange(n):
      ans = ''.join([self.rules[g] for g in ans])
    return ans
  

#this returns as much of the word (end) as we can
#conclude from what is given (the image of the cyl set) 
def sig_word_under_aut(w, A, rank=2, n=1):
  if n == 0:
    return w
  gens = alphabet[:rank] + inverse(alphabet[:rank])
  last_letter_inverse = w[-1].swapcase()
  w_appended = [w + g + h for g in gens for h in gens \
                 if g != last_letter_inverse and g!=h.swapcase()]
  A_on_w_appended = w_appended
  for j in xrange(n):
    A_on_w_appended = A.ap(A_on_w_appended)
  #print A_on_w_appended       
  j=0
  maxJ = min([len(x) for x in A_on_w_appended])
  while True:
    jsOK = True
    for k in xrange(1,len(A_on_w_appended)-1):
      #print (j, maxJ, k, A_on_w_appended)        
      if j >= maxJ or A_on_w_appended[0][j] != A_on_w_appended[k][j]:
        jsOK = False
        break
    if jsOK:
      j += 1
    else:
      break
  #print A_on_w_appended       
  return A_on_w_appended[0][:j]


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


#gives an orientation to a tripod, where the graph determining it 
# is the square, with each vertex labeled g, G in order on a loop
# each loop is length 1, for whatever reason
def orient_t(t):  
  gens = alphabet[:4]
  gensi = inverse(alphabet[:4])
  g_and_i = [g for g in gens if g in t and inverse(g) in t]
  if len(g_and_i) > 0 :
    i=t.index(g_and_i[0])
    j=t.index(inverse(g_and_i[0]))
    k = [m for m in xrange(3) if m != i and m != j][0]
    if i < k:
      i += 3
    if j < k:
      j += 3
    if i < j:
      return -1
    else:
      return 1
  else:
    tl = [x.lower() for x in t]
    if cyclic_order_coefficient('abcd', tl[0], tl[1], tl[2]) > 0:
      return -4
    else:
      return 4















