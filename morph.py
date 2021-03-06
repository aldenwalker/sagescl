import random as RAND
import fractions
import itertools

import word

from sage.all import *

alphabet = list('abcdefghijklmnopqrstuvwxyz')

class morph:
  def __init__(self, r):
    if type(r) == dict:
      self.rules = r
      newRules = {}
      for x in self.rules:
        if x.swapcase() not in self.rules:
          newRules[x.swapcase()] = word.inverse(self.rules[x])
      self.rules.update(newRules)
    else:
      self.rules = {}
      pairs = r.split(',')
      for p in pairs:
        if len(p.replace(' ','')) == 0:
          continue
        p1, p2 = p.split('->')
        self.rules[ p1.replace(' ','') ] = p2.replace(' ','')
      newRules = {}
      for x in self.rules:
        if x.swapcase() not in self.rules:
          newRules[x.swapcase()] = word.inverse(self.rules[x])
      self.rules.update(newRules)      
        
  def ap(self, L, n=1, marked=False):
    if isinstance(L, list):
      return [self.ap(x,n,marked) for x in L]
    elif L=='' or n==0:
      return ''
    else:
      ans = L
      for i in xrange(n):
        lets = list(ans)
        if marked:
          ans = word.multiply_words_marked([(self.rules[let] if let != '.' else '.') for let in lets])
        else :
          ans = word.multiply_words([self.rules[let] for let in lets])
      return ans
      
  def iterate(self, n):
    newRules = dict([ (x,x) for x in self.rules])
    for i in xrange(n):
      for x in newRules:
        newRules[x] = self.ap(newRules[x])
    return morph(newRules)
    
  def __mul__(self, other):
    if not isinstance(other, morph):
      return None
    return morph( dict( [(x,self.ap(other.ap(x))) for x in self.rules] ))
  
  def __str__(self):
    ans = ''
    for x in sorted(self.rules.keys()):
      if x.islower():
        ans += str(x) + ' -> ' + str(self.rules[x]) + ',   '
    return ans
  
  def __repr__(self):
    return str(self)

  def __call__(self, x):
    return self.ap(x)

  def __eq__(self, other):
    #print "Called with ", self, " and ", other 
    sk = set([x for x in self.rules.keys() if x.islower()])
    ok = set([x for x in other.rules.keys() if x.islower()])
    if sk != ok:
      return False
    else:
      return not any([self.rules[l] != other.rules[l] for l in sk])
  
  def is_identity(self):
    return all([self.rules[g]==g for g in self.rules])
  
  def is_conjugation(self):
    redrules = word.simultaneously_conjugated([self.rules[g] for g in self.rules if g.islower()])
    if redrules[1] == [g for g in self.rules if g.islower()]:
      return True
    else:
      return False
  
  def act_on_tripod(self, T):
    it = self.ap(T)
    min_len = min(map(len, it))
    i = 0
    while i < min_len and it[0][i] == it[1][i] and it[1][i] == it[2][i]:
      i += 1
    it = [w[i:] for w in it]
    #ensure that if any two things agree, it's 0 and 1
    if it[1][0] == it[2][0]:
      it = [it[1], it[2], it[0]]
    elif it[0][0] == it[2][0]:
      it = [it[2], it[0], it[1]]
    min_len = min_len - i
    i=0
    while i < min_len and it[0][i] == it[1][i]:
      i += 1
    return [it[0][i:], it[1][i:], word.inverse(it[0][:i]) + it[2]]
    
  def homology_matrix(self):
    rank = len(self.rules)/2
    gens = sorted([x for x in self.rules.keys() if x.islower()])
    rows = [ [ (self.rules[x].count(y) - self.rules[x].count(word.inverse(y)))  for x in gens] for y in gens]
    return rows
    
  def fixed_space(self):
    H = self.homology_matrix()
    for i in xrange(len(H)):
      H[i][i] -= 1
    #SAGE IS REQUIRED FOR THIS STUFF
    M = Matrix(H)
    V = M.right_kernel(basis='LLL')
    return list(V.basis())
    
  def action_matrix_transpose(self, n, W_in=None):
    return Matrix(self.action_matrix(n, W_in)).transpose()
  
  #returns the matrix giving the action on words of length n
  # the endomorphism must have locally blocked cancellation
  # the matrix ACTS ON THE RIGHT (ON ROW VECTORS)  
  # this does the algorithm from the traintracks paper
  def action_matrix(self, n, W_in=None):
    rank = len(self.rules)/2
    if W_in == None:
      W = all_words_of_len(n, alphabet[:rank])
    else:
      W = W_in
    M = []
    for w in W:
      targets = self.ap([w[0], w[1:]])
      Aw = self.ap(w)
      cancel_index = 0
      while targets[0][-cancel_index-1] == targets[1][cancel_index].swapcase():
        cancel_index += 1
      letters_starting_in_0 = len(targets[0]) - cancel_index 
      words_to_cancel = [ targets[1][i:i+n] for i in xrange(cancel_index) ]
      words_to_add = [ Aw[i:i+n] for i in xrange(letters_starting_in_0) ]
      M.append( [0 for wp in W] )
      for added_word in words_to_add:
        M[-1][ W.index(added_word) ] += 1
      for subbed_word in words_to_cancel:
        M[-1][ W.index(subbed_word) ] -= 1
    return M
  
      
    
  #this returns the matrix which is (I + NM), where I is the identity, 
  #N is the matrix taking a weight to its inverse, and M is the matrix which is the 
  #action of self.  Here we want a left action (on column vectors), so M will be the 
  #transpose of the output of action_matrix
  def id_plus_negative_phi(self, n, W_in=None):
    rank = len(self.rules)/2
    if W_in == None:
      W = all_words_of_len(n, alphabet[:rank])
    else:
      W = W_in
    dim = len(W)
    I = identity_matrix(dim)
    #construct the inversion matrix
    N = Matrix(dict([((W.index(inverse(W[i])), i), 1) for i in xrange(dim)]))
    #get the action matrix
    M = Matrix(self.action_matrix(n, W))
    M = M.transpose()
    return I + (N*M)    
    
  def negative_phi(self, n, W_in=None) :
    rank = len(self.rules)/2
    if W_in == None:
      W = all_words_of_len(n, alphabet[:rank])
    else:
      W = W_in
    dim = len(W)
    #construct the inversion matrix
    N = Matrix(dict([((W.index(inverse(W[i])), i), 1) for i in xrange(dim)]))
    #get the action matrix
    M = Matrix(self.action_matrix(n, W))
    M = M.transpose()
    return N*M
    
    
  def is_expanding(self):
    for g in self.rules:
      if g.isupper():
        continue
      w = self.rules[g]
      Lw = len(w)
      if Lw == 0:
        return False
      max_left_cancel = 0
      for g2 in self.rules:
        if g2.swapcase() == g:
          continue
        w2 = self.rules[g2]
        i=0
        Lw2 = len(w2)
        if Lw2 == 0:
          return False
        while w2[-(i+1)].swapcase() == w[i]:
          i += 1
          if i >= Lw:
            return False
          if i >= Lw2:
            break
        if i > max_left_cancel:
          max_left_cancel = i
      max_right_cancel = 0
      for g2 in self.rules:
        if g2.swapcase() == g:
          continue
        w2 = self.rules[g2]
        Lw2 = len(w2)
        if Lw2 == 0:
          return False
        i=0
        while w[-(i+1)] == w2[i].swapcase():
          i += 1
          if i >= Lw:
            return False
          if i >= Lw2:
            break
        if i > max_right_cancel:
          max_right_cancel = i      
      if max_left_cancel > len(w)-max_right_cancel-1:
        return False
    return True


  
def are_symmetric(A,B):
  switch = morph({'a':'b','b':'a'})
  inva = morph({'a':'a','b':'b'})
  invb = morph({'a':'a','b':'b'})
  SB = switch*B*switch
  if A == SB:
    return True
  if A == inva*B*invb or A == inva*SB*inva:
    return True
  if A == invb*B*invb or A == invb*SB*invb:
    return True
  if A == inva*invb*B*invb*inva or A == inva*invb*SB*invb*inva:
    return True
  
def power_is_reducible(A, n=2):
  ng = len(A.rules)/2
  gens = alphabet[:ng]
  for i in xrange(n):
    B = A.iterate(i+1)
    if sorted([power_reduce(cyc_red(B.ap(g)).lower())[0] for g in gens]) == alphabet[:ng]:
      return True
  return False
  
def reduce_aut_list(L):
  classes = []
  for A in L:
    found = False
    for C in classes:
      if are_symmetric(A, C):
        found = True
    if not found:
      classes.append(A)
  return classes

def random_automorphism(rank, n=None):
  gens = alphabet[:rank]
  aut = {}
  for x in gens:
    aut[x] = x
    aut[x.swapcase()] = x.swapcase()
  K = aut.keys()
  for i in xrange( (n if n!=None else 4*rank)):
    if random.random() < 0.5:  # swap
      j = random.choice(K)
      k = random.choice(K)
      tempj = aut[j]
      tempJ = aut[j.swapcase()]
      aut[j] = aut[k]
      aut[j.swapcase()] = aut[k.swapcase()]
      aut[k] = tempj
      aut[k.swapcase()] = tempJ
    else:  # multiply
      j = random.choice(K)
      k = random.choice(K)
      while k == j or k == j.swapcase():
        k = random.choice(K)
      aut[j] = multiply_words(aut[j], aut[k])
      aut[j.swapcase()] = inverse(aut[j])
  return morph(aut)


def random_composition(autos, n):
  if n==0:
    return None
  m = RAND.choice(autos)
  #print m
  for i in xrange(n-1):
    m2 = RAND.choice(autos)
    #print m2
    m = m*m2
  return m
  
def all_compositions(autos, n, ball=False):
  if ball:
    for A in itertools.chain(*[all_compositions(autos, i) for i in xrange(1,n+1)]):
      yield A
    return
  #print "Yielding length ",n
  T = Tuples(range(len(autos)), n)
  for t in T:
    ans = autos[t[0]]
    for i in xrange(1,n):
      ans = ans*autos[t[i]]
    yield ans


def random_automorphism_WH(rank, WH_gens, n=None) :
  gens = alphabet[:rank]
  aut = {}
  LWH = len(WH_gens)
  for x in gens:
    aut[x] = x
    aut[x.swapcase()] = x.swapcase()
  aut = morph(aut)
  for i in xrange( (4*rank if n==None else n) ):
    ind = random.randint(0, LWH-1)
    aut = aut*WH_gens[ind]
  return aut
  

def fixed_vector(A):
  if len(A.rules.keys()) == 4:
    H = homology_matrix(A)
    if (H[0][0]-1)*(H[1][1]-1) - H[0][1]*H[1][0] == 0:
      if [H[0][1], 1 - H[0][0] ] == [0,0]:
        if H == [[1,0],[0,1]]:
          return [1,0]
        else:
          return reduce_int_vector([H[1][1]-1, -H[1][0]])
      else:
        return reduce_int_vector([H[0][1], 1 - H[0][0] ]) 
    else:
      return None
  else:
    ### THIS WILL ONLY WORK FROM WITHIN SAGE
    M = Matrix(ZZ, homology_matrix(A))
    e_vals = M.eigenvalues()
    if 1 in e_vals:
      ev = M.eigenvectors_right()
      for e in ev:
        if e[0] == 1:
          l = lcm([x.denominator() for x in e[1][0] if x != 0])
          return l * e[1][0]  ## return just one of the eigenvectors
    return None

def random_automorphism_with_fixed_vector(rank, WH_gens, n=None) :
  while True:
    A = random_automorphism_WH(rank, WH_gens, n)
    #A = random_automorphism(rank, n)
    H = homology_matrix(A)
    if (H[0][0]-1)*(H[1][1]-1) - H[0][1]*H[1][0] == 0:
      if [H[0][1], 1 - H[0][0] ] == [0,0]:
        if H == [[1,0],[0,1]]:
          return [A, [1,0]]
        else:
          return [A, [H[1][1]-1, -H[1][0]] ]
      else:
        return [A, [H[0][1], 1 - H[0][0] ]  ]
  
def random_homomorphism_with_fixed_vector(rank, n, m) :
  while True:
    A = morph({'a':random_reduced_word(n), 'b':random_reduced_word(m)})
    H = homology_matrix(A)
    if (H[0][0]-1)*(H[1][1]-1) - H[0][1]*H[1][0] == 0:
      if [H[0][1], 1 - H[0][0] ] == [0,0]:
        if H == [[1,0],[0,1]]:
          return [A, [1,0]]
        else:
          return [A, reduce_int_vector([H[1][1]-1, -H[1][0]] ) ]
      else:
        return [A, reduce_int_vector([H[0][1], 1 - H[0][0] ]) ]


def random_homomorphism(rank, target_lens):
  targets = [random_reduced_word(n, rank) for n in target_lens]
  A= morph( dict( [ (alphabet[i], targets[i]) for i in xrange(rank)] ) )
  return A


def reduce_int_vector(L):
  g = reduce(fractions.gcd, L)
  return [x/g for x in L]


def homology_matrix(A):
  rank = len(A.rules)/2
  gens = sorted([x for x in A.rules.keys() if x.islower()])
  rows = [ [ (A.rules[x].count(y) - A.rules[x].count(inverse(y)))  for x in gens] for y in gens]
  return rows

class whiteheadAut(morph):
  def __init__(self, rank, A, a):
    self.A = A
    self.a = a
    self.rules = {}
    for gen in alphabet[:rank] + map(word.inverse, alphabet[:rank]):
      if gen == a or gen == word.inverse(a):
        self.rules[gen] = gen
      elif gen in A and word.inverse(gen) not in A:
        self.rules[gen] = gen+a
      elif gen not in A and word.inverse(gen) in A:
        self.rules[gen] = word.inverse(a)+gen
      elif gen not in A and word.inverse(gen) not in A:
        self.rules[gen] = gen
      else:
        self.rules[gen] =word.inverse(a) + gen + a


def WAgen(rank):
  if rank == 0:
    return
  allGens = alphabet[:rank] + map(word.inverse, alphabet[:rank])
  LA = len(allGens)
  #first give the identity, then never do it again
  yield whiteheadAut(rank, [], allGens[0])
  
  #iterate over the multiplier
  for a in allGens:
    #the number of gens possibilities to go into A is 2^(len(allGens)-2)
    for Anum in xrange(1,2**(LA-2)): #0 would be the identity
      currentAIndex = 0
      currentGenIndex = 0
      A = []
      while currentGenIndex < LA:
        if a.lower() == allGens[currentGenIndex].lower():
          currentGenIndex += 1
          continue
        if (Anum >> currentAIndex)&1 == 1:
          A.append(allGens[currentGenIndex])
        currentGenIndex += 1
        currentAIndex += 1
      yield whiteheadAut(rank, A, a)

def IAgen(rank):
  if rank < 3:
    return
  gens = alphabet[:rank]
  semigroup_gens = alphabet[:rank] + map(word.inverse, alphabet[:rank])
  #inverses:
  #for i in xrange(rank):
    #d = dict([(ell, (ell if j!=i else word.inverse(ell))) for j,ell in enumerate(alphabet[:rank])])
    #yield morph(d)
  #permutations:
  #for i in xrange(rank):
    #for j in xrange(rank):
      #if j==i:
        #continue
      #d = dict([(ell, (alphabet[i] if k==j else (alphabet[j] if k==i else ell))) for k,ell in enumerate(alphabet[:rank])])
      #yield morph(d)
  #conjugations
  for i in xrange(rank):
    g = alphabet[i]
    for j in xrange(rank):
      if j == i:
        continue
      d = dict([(ell, (g+ell+word.inverse(g) if k==j else ell)) for k,ell in enumerate(alphabet[:rank])])
      yield morph(d)
  #multiply by commutator of other gens
  #throw in inverse commutators (in the other order) too
  for i in xrange(rank):
    g1 = alphabet[i]
    for j in xrange(rank):
      if j==i:
        continue
      g2 = alphabet[j]
      c = g1+g2+word.inverse(g1)+word.inverse(g2)
      for k in xrange(rank):
        if k == i or k == j:
          continue
        d = dict([(ell, (ell+c if m==k else ell)) for m,ell in enumerate(alphabet[:rank])])
        yield morph(d)
  

def random_IA(rank, n, IAg=None):
  if IAg == None:
    IAg = list(IAgen(rank))
  ans = morph(dict([(ell,ell) for ell in alphabet[:rank]]))
  LIAg = len(IAg)
  for i in xrange(n):
    gi = RAND.randint(0,LIAg-1)
    ans = ans*IAg[gi]
  return ans




# assumes rank of whatever it finds
def min_in_orbit(c, R=None, m=False):
  #make a copy of c
  if isinstance(c, str):
    C = [c]
  else:
    C = [x for x in c]
    
  if R == None:
    rank = word.chain_rank(c)
  else:
    rank = R
  
  Clen = sum(map(len, C))
  
  #go through all whitehead auts, restarting if we get better
  while True:
    WA = WAgen(rank)
    gotBetter = False
    for W in WA:
      newC = W.ap(C, marked=m)
      newLen = sum(map(len,newC))
      if newLen < Clen:
        C = newC
        Clen = newLen
        gotBetter = True
        break
    if not gotBetter:
      break
  
  if isinstance(c, str):
    return C[0]
  else:
    return C

def is_minimal_in_orbit(c, R=None, m=False):
  if isinstance(c, str):
    C = [c]
  else:
    C = c
  if R == None:
    rank = word.chain_rank(C)
  else:
    rank = R
  CLen = sum(map(len, C))
  WA = WAgen(rank)
  for W in WA:
    newC = W.ap(C, marked=m)
    newLen = sum(map(len, newC))
    if newLen < CLen:
      return false
  return true

def get_minimal_elements(L, R=None, m=False):
  M = []
  WA = list(WAgen(R))
  for w in L:
    wL = len(w)
    minimal = True
    for W in WA:
      Ww = W.ap(w,marked=m)
      if len(Ww) < wL:
        minimal = False
        break
    if minimal:
      M.append(w)
  return M
