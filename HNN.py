#!/usr/bin/python

import fractions
import subprocess
import random as RAND
import sys
import os

from morph import *
#load('morph.py')
from scl import *
from word import *
from sage.all import *

import covering
import fatgraph

alphabet = list('abcdefghijklmnopqrstuvwxyz')


def lcm(a,b):
  return (a*b)/fractions.gcd(a,b)
  
def integerize_fractions(L):
  denoms = [x.denominator for x in L]
  g = reduce(lcm, denoms)
  cleared = [int(g*x) for x in L]
  g = reduce(fractions.gcd, cleared)
  return [x/g for x in cleared]
  
def immersed_surface(CO, C1_in):
  """Gives an immersed surface bounding bd(S) + C1, where bd(S) is the 
  boundary of the realization with cyclic order CO"""
  C1 = (C1_in if type(C1_in)==list else [C1_in])
  bd = boundary(CO)
  pt = smallest_positive_scl_vertex(bd, C1)
  full_chain = bd + C1
  weights = [pt[0] for x in bd] + [pt[1] for x in C1]
  #simplify the weights
  weights = integerize_fractions(weights)
  return extremal_surface(full_chain, weights)
  
def realization_boundary_patterns_of_immersed_surface(CO, C1_in):
  surf = immersed_surface(CO, C1_in)
  cyclic_boundary = boundary(CO)
  all_patterns = surf.boundaries_with_patterns()
  return [x for x in all_patterns if \
          cyclically_contained(power_reduce(x)[0], cyclic_boundary)]
  



def build_silly_surface(w):
  #build the surface (a stupid surface)
  assigned_letters = [[0 for x in W] for W in w]
  current_junctions = []
  lw = len(w)
  all_letters = [(i,j) for i in xrange(lw) for j in xrange(len(w[i]))]
  #make the pairings
  pairs = {}
  for (i,j) in all_letters:
    if assigned_letters[i][j] == 1:
      continue
    assigned_letters[i][j] = 1
    for (k,l) in all_letters:
      if assigned_letters[k][l] == 0 and w[i][j] == w[k][l].swapcase():
        assigned_letters[k][l] = 1
        pairs[(i,j)] = (k,l)
        pairs[(k,l)] = (i,j) 
        break
  #print pairs
  #figure out the sequences
  while len(all_letters) > 0:
    start = all_letters[0]
    current_letter = start
    all_letters.remove(start)
    current_junctions.append([ (start, pairs[start]) ])
    #print current_junctions
    #print all_letters 
    while True:
      paired_letter = pairs[current_letter]
      (i,j) = paired_letter
      next_letter = (i, (j+1)%len(w[i]))
      if next_letter == start:
        break
      current_letter = next_letter
      current_junctions[-1].append( (current_letter, pairs[current_letter]) )
      all_letters.remove(current_letter)
  return current_junctions

#note the pairs go backwards because we need to 
#reverse the word (it's a left action)
#also, we need to antisymmetrize it
def tripod_rot(CO, t):
  p1 = pair_jump(t[1], inverse(t[0]), CO) - pair_jump(t[0], inverse(t[1]), CO)
  p2 = pair_jump(t[2], inverse(t[1]), CO) - pair_jump(t[1], inverse(t[2]), CO)
  p3 = pair_jump(t[0], inverse(t[2]), CO) - pair_jump(t[2], inverse(t[0]), CO)
  return p1 + p2 + p3

def rot_power(CO, w_in, TA, n):
  if type(w_in) == str:
    w = [cyc_red(w_in)]
  else:
    w = cyc_red(w_in)
  rank = chain_rank(w)
  #figure out what to add to make it homologically trivial
  for g in alphabet[:rank]:
    c = sum([x.count(g) - x.count(g.upper()) for x in w])
    if c > 0:
      w.append(c*g.upper())
    elif c < 0:
      w.append((-c)*g)
  #get the surface
  S = build_silly_surface(w)
  #print S
  #get the tripods
  tri_list = []
  for s in S:
    for i in xrange(1,len(s)-1):
      (i1,j1) = s[0][0]
      (i2,j2) = s[i][0]
      (i3,j3) = s[i+1][0]
      tri_list.append( (w[i1][j1], w[i2][j2], w[i3][j3]) )
  #print tri_list
  tri_list = [t for t in tri_list if len(set(t))==3]
  #print tri_list
  #act on the tripods
  tri_list = TA.ap(tri_list, n)
  #get rot
  r = sum([tripod_rot(CO, t) for t in tri_list])
  if r%(4*rank) != 0:
    print "not divisible by 4*rank?"
  return r / (4*rank)
  



def boundary(w):
  wl = len(w)
  boundaries = []
  found = [0 for i in xrange(wl)]
  while True:
    try:
      ind = found.index(0)
    except:
      break 
    let = w[ind]
    letI = inverse(let)
    found[ind] = 1
    indI = w.index(letI)
    boundaries.append(let)
    ind = (indI+1)%wl
    while found[ind] != 1:
      let = w[ind]
      letI = inverse(let)
      indI = w.index(letI)
      found[ind] = 1
      boundaries[-1] += let
      ind = (indI+1)%wl
      
  return boundaries
      
    


def hillclimb_orders(rank, f, restarts, iters):
  gens = alphabet[:rank] + inverse(alphabet[:rank])
  best_found = None
  best_score = 0
  L = 2*rank
  for i in xrange(restarts):
    current_order = [x for x in gens]
    potential_step = [x for x in gens]
    random.shuffle(current_order)
    score = f(''.join(current_order))
    #print "restart"; sys.stdout.flush()
    for j in xrange(iters):
      broke = False
      for s1 in xrange(L):
        for s2 in xrange(s1+1, L):
      #s1 = random.randint(0,L-1)
      #s2 = random.randint(0,L-1)
      #while s2 == s1:
      #  s2 = random.randint(0,L-1)
          temp = potential_step[s1]
          potential_step[s1] = potential_step[s2]
          potential_step[s2] = temp
          potential_score = f(''.join(potential_step))
          if potential_score > score:
            current_order = [x for x in potential_step]
            score = potential_score
            broke = True
            if score > best_score:
              best_score = score
              best_found = [x for x in current_order]
              print [best_score, ''.join(best_found)]
              sys.stdout.flush()
            break
          else:
            temp = potential_step[s1]
            potential_step[s1] = potential_step[s2]
            potential_step[s2] = temp
        if broke:
          break
      if not broke:
        break
  return [best_score, ''.join(best_found)]
      


def random_cyclic_order(rank):
  gens = alphabet[:rank] + inverse(alphabet[:rank])
  random.shuffle(gens)
  return ''.join(gens)
  

 
def word_from_vector(v):
  lets =  [ (v[k]*alphabet[k] if v[k] >=0 \
                       else (-v[k])*inverse(alphabet[k]))      \
                            for k in xrange(len(v)) ]
  return ''.join(lets)


def chain_from_vector(v):
  lets =  [ (v[k]*alphabet[k] if v[k] >=0 \
                       else (-v[k])*inverse(alphabet[k]))      \
                            for k in xrange(len(v)) ]
  return lets

  
#is the triple w1 allowed to be paired with w2?
def triple_pairing_allowed(w1, w2, CO):
  if w1[1] != w2[1].swapcase():
    return False
  t1 = (w1[0].swapcase(), w1[1], w2[2])
  t2 = (w1[2], w2[0].swapcase(), w2[1])
  if t1[0] == t1[1] or t1[0] == t1[2] or t1[1] == t1[2]:
    pass 
  elif tripod_sign(CO, t1) == -1:
    return False
  if t2[0] == t2[1] or t2[0] == t2[2] or t2[1] == t2[2]:
    pass 
  elif tripod_sign(CO, t2) == -1:
    return False
  return True


def generate_triple_pairs(CO):
  rank = len(CO)/2
  W = all_words_of_len(3, alphabet[:rank])
  LW = len(W)
  pairs = []
  for i in xrange(LW):
    for j in xrange(i+1,LW):
      if triple_pairing_allowed(W[i], W[j], CO):
        pairs.append( (W[i], W[j]) )
  return pairs

#generate rows which must all be orthogonal to the vector of words
#of length 3, plus the pairs
#the pairs are in order
def generate_pair_constraints(CO, pairs):
  rank = len(CO)/2
  W = all_words_of_len(3, alphabet[:rank])
  LW = len(W)
  LP = len(pairs)
  rows = []
  dim = LW + LP
  for i in xrange(LW):
    this_row = dim * [0]
    this_row[i] = -1
    for j in xrange(LP):
      if W[i] in pairs[j]:
        this_row[LW + j] = 1
    rows.append(this_row)
  return rows
    

def generate_w_ell_cons(rank, ell, pairs):
  verts = all_words_of_len(ell-1, alphabet[:rank])
  W = all_words_of_len(ell, alphabet[:rank])
  LW = len(W)
  LP = len(pairs)
  dim = LW + LP
  rows = []
  for i in xrange(len(verts)):
    this_row = dim * [0]
    for j in xrange(LW):
      if W[j][1:] == verts[i]:
        this_row[j] += -1
      if W[j][:2] == verts[i]:
        this_row[j] += 1
    rows.append(this_row)
  return rows

def inverse_matrix(rank, ell):
  W = all_words_of_len(ell, alphabet[:rank])
  LW = len(W)
  rows = []
  for i in xrange(LW):
    this_row = LW*[0]
    this_row[ W.index(inverse(W[i])) ] = 1
    rows.append(this_row)
  return rows
  
def next_permutation(p_in):
  p = [x for x in p_in]
  pl = len(p)
  i = pl-2
  while i>=0 and p[i]>p[i+1]:
    i -= 1
  if i == -1:
    return None
  j=i+1
  while j<pl and p[j]>p[i]:
    j += 1
  j -= 1
  temp = p[i]
  p[i] = p[j]
  p[j] = temp
  p[i+1:] = p[i+1:][::-1]   #reverse after i
  return p
  
  

def permutations(limit):
  p = [i for i in xrange(limit)]
  P = [p]
  while True:
    p = next_permutation(p)
    if p == None:
      break 
    P.append(p)
  return P



def next_cyclic_order(CO):
  rank = len(CO)/2
  gens = alphabet[:rank] + inverse(alphabet[:rank])
  indices = dict([ (gens[i], i) for i in xrange(len(gens)) ])
  P = [indices[g] for g in CO[1:]]
  #print P,
  P = next_permutation(P)
  #print P
  if P == None:
    return None
  return ''.join([CO[0]] + [gens[p] for p in P])

def cyclic_orders(rank):
  gens = alphabet[:rank]
  gens += inverse(gens)
  P = permutations(2*rank-1)
  CO = [ [gens[0]] + [gens[p[j]+1] for j in xrange(2*rank-1)] for p in P]
  return [''.join(x) for x in CO]


def all_tripods(rank):
  letters = alphabet[:rank]
  letters += inverse(letters)
  AT = []
  for i in xrange(2*rank):
    l0 = letters[i]
    AT.extend( [ (l0, letters[j], letters[k]) for j in xrange(i+1, 2*rank) \
                                              for k in xrange(i+1, 2*rank) if j!= k ] )
  return AT
  
  

def simplify_tripod(T):
  if any([len(x)==0 for x in T]):
    return None
  l1,l2,l3 = [x[0] for x in T]
  if l1 == l2:
    if l1 == l3:
      return simplify_tripod([x[1:] for x in T])
    else:
      return simplify_tripod([ T[0][1:], T[1][1:], inverse(T[0][0]) + T[2]])
  if l1 == l3:
    return simplify_tripod([ T[0][1:], inverse(T[0][0]) + T[1], T[2][1:] ] )
  if l2 == l3:
    return simplify_tripod([ inverse(T[1][0]) + T[0], T[1][1:], T[2][1:] ] )
  return tuple([x[0] for x in T])

 
#this is only rank 2
def recognize_tripod(T):
  which = [x for x in ['a','b','A','B'] if x not in T][0]
  ind_of_inverse = T.index(inverse(which))
  rotated = [T[(ind_of_inverse+i)%3] for i in xrange(3)]
  if rotated[1].isupper():
    return (-1, which)
  else:
    return (1,which)


  
def tripod_canonical_form(T):
  swapped = [x.swapcase() for x in T]
  m = min(swapped)
  i = swapped.index(m)
  return T[i:] + T[:i]
  
  

def tripod_flip(T):
  return (T[0], T[1], T[2])

def tripod_sign(CO, T):
  'Is a tripod cyclically oriented with respect to the cyclic order CO?'
  inds = [CO.index(t) for t in T]
  if inds[1] < inds[0]:
    inds[1] += len(CO)
  if inds[2] < inds[0]:
    inds[2] += len(CO)
  return (1 if inds[1] < inds[2] else -1)

def order_to_tripods(CO):
  'Return an exhausting set of tripods'
  return [ (CO[0], CO[i], CO[i+1]) for i in xrange(1, len(CO)-1)]

def act_on_tripod(A, tripod):
  if type(tripod) ==list:
    return [A.act_on_tripod(x) for x in tripod] 
  elif type(tripod) == tuple:
    T = tripod
  #print actual_tripod
  acted_on = tuple(A.ap(list(T)))
  result = simplify_tripod(acted_on)
  return result

def tripod_action(A):
  T = all_tripods(len(A.rules)/2)
  act = [ (t, act_on_tripod(A,t)) for t in T]
  try:
    return tripod_morph(dict(act))
  except TypeError:
    return None

class tripod_morph:
  def __init__(self, r):
    self.rules = dict([ (tripod_canonical_form(x), tripod_canonical_form(r[x])) for x in r])
    newRules = {}
    for x in self.rules:
      if tripod_flip(x) not in self.rules:
        newRules[tripod_flip(x)] = tripod_flip(self.rules[x])
    self.rules.update(newRules)
  
  def ap(self, T, n=1):
    if type(T) == list:
      return [self.ap(t, n) for t in T]
    else:
      tripod = tripod_canonical_form(T)
    
    for i in xrange(n):
      tripod = self.rules[tripod]
    
    return tripod
  
  def preserves_order(self, CO, n=1):
    T = order_to_tripods(CO)
    return all([ tripod_sign(CO, t)==1 for t in self.ap(T, n)])
    

def to_mathematica(L):
  return str(L).replace('[','{').replace(']','}')
      

class whiteheadAut(morph):
  def __init__(self, rank, A, a):
    self.A = A
    self.a = a
    self.rules = {}
    for gen in alphabet[:rank] + map(inverse, alphabet[:rank]):
      if gen == a or gen == inverse(a):
        self.rules[gen] = gen
      elif gen in A and inverse(gen) not in A:
        self.rules[gen] = gen+a
      elif gen not in A and inverse(gen) in A:
        self.rules[gen] = inverse(a)+gen
      elif gen not in A and inverse(gen) not in A:
        self.rules[gen] = gen
      else:
        self.rules[gen] = inverse(a) + gen + a


def WAgen(rank):
  allGens = alphabet[:rank] + map(inverse, alphabet[:rank])
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

# assumes rank of whatever it finds
def min_in_orbit(c, R=None):
  #make a copy of c
  if isinstance(c, str):
    C = [c]
  else:
    C = [x for x in c]
    
  if R == None:
    rank = 0
    for w in C:
      for let in w:
        ind = alphabet.index(let.lower())
        if ind+1 > rank:
          rank = ind+1
  else:
    rank = R
  
  Clen = sum(map(len, C))
  
  #go through all whitehead auts, restarting if we get better
  while True:
    WA = WAgen(rank)
    gotBetter = False
    for W in WA:
      newC = W.ap(C)
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

def gen_tuples(limits):
  return reduce( lambda current, next_limit:                                 \
                    [x + [i] for x in current for i in xrange(next_limit)],  \
                    limits,                                                  \
                    [[]] )

#this modifies the input!
def next_tuple(current, limits):
  lc = len(current)
  i = lc-1
  next_t = [x for x in current]
  while i >= 0 and current[i] == limits[i]-1:
    i -= 1
  if i == -1:
    return None
  next_t[i] += 1
  for j in xrange(i+1, lc):
    next_t[j] = 0
  return next_t
    
def tuples_gen(limits):
  ll = len(limits)
  t = [0 for i in xrange(ll)]
  while t != None:
    yield t
    t = next_tuple(t, limits)
  

def stupid_endomorphism(targets):
  if any([len(w)==0 for w in targets]):
    return True
  
  pairs = [(i,j) for i in xrange(len(targets)) for j in xrange(len(targets)) \
                 if i!=j]
  if any([targets[i]==targets[j] for (i,j) in pairs]):
    return True
  
  if sorted([cyc_red(w).lower() for w in targets]) == alphabet[:len(targets)]:
    return True
    
  if any([cyc_red(targets[k]).lower() == alphabet[k] for k in xrange(len(targets))]):
    return True

  return False
  
  

  
#this writes a matrix to a file, in the format expected by trollop
def write_matrix(f, M, words_in_order):
  f.write(str(len(M)) + ' ' + str(len(M[0])) + '\n')
  f.write(' '.join(words_in_order) + '\n')
  for row in M:
    f.write(' '.join([str(x) for x in row]) + '\n')

def write_vector(f, v):
  f.write(str(len(v)) + '\n')
  f.write(' '.join([str(x) for x in v]) + '\n')

#this writes out the two matrices and vector, as expected by the matrix function
#of trollop
def traintrack_matrices(directory, A, n, v=None, filename=''):
  rank = len(A.rules)/2
  words = all_words_of_len(n, alphabet[:rank])
  IdpNP = list(A.id_plus_negative_phi(n, words))
  if v == None:
    v = list(A.fixed_space())[0]
  hom_matrix = [ [(sign(x[1]) if x[1].lower() == g else 0) for x in words] for g in alphabet[:rank] ]
  Mfile = directory + '/' + filename + 'M.mat'
  Nfile = directory + '/' + filename + 'N.mat'
  vfile = directory + '/' + filename + 'b.vec'
  fM = open(Mfile,'w')
  fN = open(Nfile,'w')
  fv = open(vfile,'w')
  write_matrix(fM, IdpNP, words)
  write_matrix(fN, hom_matrix, words)
  write_vector(fv, v)
  fM.close()
  fN.close()
  fv.close()
  
#this writes out the two matrices and vector, as expected by the matrix function
#of trollop
def traintrack_matrices_sep_domain(directory, A, n, v=None, filename=''):
  rank = len(A.rules)/2
  words = all_words_of_len(n, alphabet[:rank])
  phi = list(A.negative_phi(n, words))
  if v == None:
    v = list(A.fixed_space())[0]
  hom_matrix = [ [(sign(x[1]) if x[1].lower() == g else 0) for x in words] for g in alphabet[:rank] ]
  Mfile = directory + '/' + filename + 'M.mat'
  Nfile = directory + '/' + filename + 'N.mat'
  vfile = directory + '/' + filename + 'b.vec'
  fM = open(Mfile,'w')
  fN = open(Nfile,'w')
  fv = open(vfile,'w')
  write_matrix(fM, phi, words)
  write_matrix(fN, hom_matrix, words)
  write_vector(fv, v)
  fM.close()
  fN.close()
  fv.close()
  
  
  
  
  
#take a word with capitals into its gap form
def to_gap(s):
  sl = list(s)
  sl = [(g if g.islower() else g.lower() + '^-1') for g in sl]
  return '*'.join(sl)
  
  
  
#use gap to find finite index subgroups
def gap_finite_index_subgroups_with_homology(A, index_bound, iteratively=False):
  rank = len(A.rules)/2
  gens = alphabet[:rank]
  relators = [to_gap(multiply_words([inverse(g), 't', A.rules[g], 'T'])) for g in gens]
  gap.eval('F:=FreeGroup(' + ','.join(['"' + g + '"' for g in gens]) + ',"t");')
  gap.eval('AssignGeneratorVariables(F);')
  gap.eval('relators:=[' + ','.join(relators) + '];')
  gap.eval('G:=F/relators;')
  betti1 = eval(gap.eval('AbelianInvariants(G);')).count(0)
  betti2 = len(A.fixed_space())
  if betti2 > 0:
    return (-1, betti2)
  chiG = 1-betti1+betti2
  #print chiG
  index_bounds_to_check = (range(4, index_bound+1) if iteratively else [index_bound])
  for i in index_bounds_to_check:
    gap.eval('subgroups:=LowIndexSubgroupsFpGroup(G,' + str(i) + ');')
    num_subgroups = int(gap.eval('Length(subgroups);'))
    subgroups = []
    for i in xrange(num_subgroups):
      index = int(gap.eval('Index(G,subgroups[' + str(i+1) + ']);'))
      ab_invar = eval(gap.eval('AbelianInvariants(subgroups[' + str(i+1) + ']);'))
      #print ab_invar
      betti1_cover = ab_invar.count(0)
      chi_cover = index * chiG
      betti2_cover = chi_cover - 1 + betti1_cover
      if betti2_cover < 0:
        print "Problem with negative betti2"
      if betti2_cover > 0:
        subgroups.append( (i, index, betti2_cover) )
    if len(subgroups) > 0:
      break
  return subgroups
  

def finite_index_subgroups_with_homology(A, indexes):
  gens = A.rules.keys()
  gens = [g for g in gens if g.islower()]
  subgroups_with_homology = []
  for ind in indexes:
    perms = permutations(ind)
    T = tuples_gen([len(perms) for i in xrange(len(gens))])
    print "Running index: ", ind
    count = 0
    total = len(perms) ** len(gens)
    for t in T:
      count += 1
      if count%1000 == 0:
        print "\r", count, " of ", total,
        sys.stdout.flush()
      selected_perms = [perms[k] for k in t]
      #print "Trying covering with perms: ", selected_perms
      G = covering.FISubgroup(gens, selected_perms)          #create the finite index subgroup
      if not G.connected:                                    #determine if it's connected
        #print "Not connected"
        continue
      G_gens_under_A = A.ap(G.gens_in_base_group())           
      if any([not G.contains(g) for g in G_gens_under_A]):   #determine if it's fixed by A
        #print "Not fixed by endo"
        continue
      AG = G.hom_from_base_hom(A)                            #get the endo on G
      #HA = AG.homology_matrix()                             #get the homology action
      V = AG.fixed_space()                                   #get the space of fixed vectors
      #print "Fixed space: ", V
      if len(V) > 0:
        subgroups_with_homology.append( (G, AG, V) )
        print subgroups_with_homology[-1]
        sys.stdout.flush()
  return subgroups_with_homology

def check_if_homology_is_geometric(A, V, CO_to_check, power_to_check=1):
  if not A.is_expanding():
    return None
  AT = tripod_action(A)
  preserved_orders = []
  immersed_surfaces
  for co in CO_to_check:
    for n in xrange(1, power_to_check+1):
      if AT.preserves_order(co):
        preserved_orders.append(co)
  for co in preserved_orders:
    #get the minimal traintrack  with th homology class
    #determine what the value of rot is on that traintrack
    #if they are equal, append
    traintrack = ''
    r=0
    if True:
      immersed_traintracks.append( (traintrack, r) )
  return immsersed_traintracks


def check_HNN_extensions_for_homology_covers(rank, word_len, index_check, only_check_expanding=False):
  W = all_words_of_len_le(word_len, alphabet[:rank])
  LW = len(W)
  T = tuples_gen([LW for i in xrange(rank)])
  all_endos_count = 0
  nonstupid_expanding_count = 0
  trivially_good_count = 0
  good_count = 0
  examples = []
  non_examples = []
  for t in T:
    targets = [W[i] for i in t]
    all_endos_count += 1
    if stupid_endomorphism(targets):
      continue
    A = morph( dict( [(alphabet[i], targets[i]) for i in xrange(rank)]) )
    if only_check_expanding and (not A.is_expanding()):
      continue
    nonstupid_expanding_count += 1
    print "Checking ", A, ' ...',
    subgroups_with_homology = gap_finite_index_subgroups_with_homology(A, index_check, iteratively=True)
    print ' done'
    sys.stdout.flush()
    if len(subgroups_with_homology) > 0:
      if subgroups_with_homology == (-1,1):
        trivially_good_count += 1
      examples.append( (A, subgroups_with_homology) )
      good_count += 1
      print "Good! - ", subgroups_with_homology[:7]
    else:
      non_examples.append( A )
  return (all_endos_count, nonstupid_expanding_count, trivially_good_count, good_count, examples, non_examples)
  
  
def check_HNN_extensions_for_ffolded(rank, word_len):
  W = all_words_of_len_le(word_len, alphabet[:rank])
  LW = len(W)
  T = tuples_gen([LW for i in xrange(rank)]) 
  ns_endo_count = 0
  ns_expanding_fixed_count = 0
  ffolded_count = 0
  found_one = False
  good_list = []
  bad_list = []
  for t in T:
    targets = [W[i] for i in t]
    if stupid_endomorphism(targets):
      continue
    ns_endo_count += 1
    A = morph( dict( [(alphabet[i], targets[i]) for i in xrange(rank)]) )
    if not A.is_expanding():
      continue
    V = A.fixed_space()
    if len(V) == 0:
      continue
    ns_expanding_fixed_count += 1
    v = V[0]
    CL = [chain_from_vector(v), [word_from_vector(v)]]
    print "Trying ", A
    found_one = False
    for C in CL:
      if is_inverse_chain(C, inverse(A.ap(C))):
        print "torus"
        continue
      CC = C + inverse(A.ap(C,marked=True), marked=True)
      print "With chain: ", C, ", ", CC
      if gallop('rose' + str(rank) + '.fg', CC, only_check_exists=True, folded=True, ffolded=True):
        print "It's good!"
        ffolded_count += 1
        good_list.append( (A,v,C) )
        found_one = True
        break
      else:
        print "No good"
        bad_list.append( (A,v,C) )
    if not found_one:
      AA = A.iterate(2)
      print "Trying the square: ", AA
      for C in CL:
        if is_inverse_chain(C, inverse(AA.ap(C))):
          print "torus"
          continue
        CC = C + inverse(AA.ap(C,marked=True), marked=True)
        if sum(map(len, CC)) > 90:
          continue
        print "With chain: ", C, ", ", CC
        if gallop('rose' + str(rank) + '.fg', CC, only_check_exists=True, folded=True, ffolded=True):
          print "It's good!"
          ffolded_count += 1
          good_list.append( (AA,v,C) )
          foudn_one = True
          break
        else:
          print "No good"
          bad_list.append( (AA,v,C) )
  return (ns_endo_count, ns_expanding_fixed_count, ffolded_count, good_list, bad_list)
        

  
def check_all_HNN_extensions_for_ffolded(rank, word_len):
  W = all_words_of_len_le(word_len, alphabet[:rank])
  LW = len(W)
  T = tuples_gen([LW for i in xrange(rank)]) 
  ns_endo_count = 0
  ffolded_count = 0
  found_one = False
  good_list = []
  bad_list = []
  for t in T:
    targets = [W[i] for i in t]
    if stupid_endomorphism(targets):
      continue
    ns_endo_count += 1
    A = morph( dict( [(alphabet[i], targets[i]) for i in xrange(rank)]) )
    CL = [random_hom_triv_word(6) for i in xrange(4)]
    print "Trying ", A
    found_one = False
    for C in CL:
      if is_inverse_chain(C, inverse(A.ap(C))):
        print "torus"
        continue
      CC = C + inverse(A.ap(C,marked=True), marked=True)
      print "With chain: ", C, ", ", CC
      if gallop('rose' + rank + '.fg', CC, only_check_exists=True, folded=True, ffolded=True):
        print "It's good!"
        ffolded_count += 1
        good_list.append( (A,v,C) )
        found_one = True
        break
      else:
        print "No good"
        bad_list.append( (A,v,C) )
  

def check_one_HNN_extension_for_ffolded(A, trials, word_len, use_fixed=True, use_hom_triv=False):
  rank = len(A.rules)/2
  fs = []
  if use_fixed:
    fs.extend([list(x) for x in A.fixed_space()])
  if use_hom_triv:
    fs.append([0 for i in xrange(rank)])
  num_vecs = len(fs)
  for i in xrange(trials):
    ind = RAND.choice(xrange(num_vecs))
    full_len = word_len + sum([abs(x) for x in fs[ind]])
    w = cyc_red(random_word_with_hom_image(full_len, rank, fs[ind]))
    C = [w, inverse(A.ap(w, marked=True), marked=True)]
    print "Trying vector ", fs[ind], " with chain: ", C
    sys.stdout.flush()
    if gallop('rose' + str(rank) + '.fg', C, only_check_exists=True, folded=True, ffolded=True):
      print "It's good!"
      return w

  
def check_chain_for_folded(C, rank=2, only_trivalent=False):
  gf = "rose" + str(rank) + ".fg"
  return gallop(gf, C, folded=True, only_check_exists=True, trivalent=only_trivalent)

def bounds_folded_stats(word_len, rank, trials=None, trivalent=False):
  if trials == None:
    W = all_words_of_len_le(word_len, alphabet[:rank])
    Wlen = len(W)
    good_count = 0
    total_count = 0
    for w in W:
      if w=='' or (not is_hom_triv(w)):
        continue
      bounds_folded = check_chain_for_folded(w, rank, trivalent)
      if total_count%100 == 0:
        print "\r", total_count, " of ", Wlen,
        sys.stdout.flush()
      if bounds_folded:
        good_count += 1
      total_count += 1
    print ""
    return (good_count, total_count)

  else:
    good_count = 0
    total_count = 0
    for i in xrange(trials):
      w = random_hom_triv_word(word_len, rank)
      bounds_folded = check_chain_for_folded(w, rank, trivalent)
      if i%100 == 0:
        print "\r", i, " of ", trials,
        sys.stdout.flush()
      if bounds_folded:
        good_count += 1
      total_count += 1
    print ""
    return (good_count, total_count)
    
    


  
def find_flat_surfaces(rank, word_len):
  W = all_words_of_len_le(word_len, alphabet[:rank])
  LW = len(W)
  T = tuples_gen([LW for i in xrange(rank)])
  all_endos_count = 0
  nonstupid_expanding_fix_hom_endos_count = 0
  nontrivial_lowerbound_count = 0
  good_examples_count = 0
  good_power_count = 0
  nontrivial_rots_count = 0
  examples = []
  CO = cyclic_orders(rank)
  for t in T:
    targets = [W[i] for i in t]
    all_endos_count += 1
    if stupid_endomorphism(targets):
      continue
    A = morph( dict( [(alphabet[i], targets[i]) for i in xrange(rank)]) )
    if not A.is_expanding():
      continue
    V = A.fixed_space()
    if len(V) == 0:
      continue
    v = V[0]
    AT = tripod_action(A)
    nonstupid_expanding_fix_hom_endos_count += 1
    good_power = 1
    good_orders = []
    while True:
      good_orders = [c for c in CO if AT.preserves_order(c, good_power)]
      if len(good_orders)>0 or good_power == 4:
        break
      else:
        good_power += 1
    if len(good_orders) == 0:
      continue 
    good_power_count += 1
    if good_power > 2:
      continue
    C = chain_from_vector(v)
    iC = cyc_red(A.ap(inverse(C), good_power))
    rots = [(c, rot(c, C+iC)) for c in good_orders if rot(c, C+iC) > 0]
    if len(rots) == 0:
      continue
    nontrivial_rots_count += 1
    s = scl(C+iC, scylla=True, scylla_i=True)
    print "Trying ", A, " on chain ", C, " with rots ", rots, " and scl=", s
    sys.stdout.flush()
    if any([2*s == x[1] for x in rots]):
      examples.append([A, C, rots])
      good_examples_count += 1
      print "Good!"
  return (all_endos_count, \
          nonstupid_expanding_fix_hom_endos_count, \
          good_power_count, \
          nontrivial_rots_count, \
          good_examples_count, \
          examples)
      
  

def find_realized_surfaces(rank, word_len, ell):
  W = all_words_of_len_le(word_len, alphabet[:rank])
  LW = len(W)
  T = tuples_gen([LW for i in xrange(rank)])
  all_endos_count = 0
  nonstupid_expanding_fix_hom_endos_count = 0
  nontrivial_lowerbound_count = 0
  good_examples_count = 0
  examples = []
  for t in T:
    targets = [W[i] for i in t]
    all_endos_count += 1
    if stupid_endomorphism(targets):
      continue
    A = morph( dict( [(alphabet[i], targets[i]) for i in xrange(rank)]) )
    if not A.is_expanding():
      continue
    V = A.fixed_space()
    if len(V) == 0:
      continue
    nonstupid_expanding_fix_hom_endos_count += 1
    for v in V:
      C = chain_from_vector(v)
      CC = C + A.ap(inverse(C))
      sclCC = scl(CC, scylla=True, scylla_i=True)
      C2 = [word_from_vector(v)]
      CC2 = C2 + A.ap(inverse(C2))
      sclCC2 = scl(CC2, scylla=True, scylla_i=True)
      #CC3 = C + A.ap(inverse(C), 2)
      #sclCC3 = scl(CC3, scylla=True, scylla_i=True)
      #CC4 = C + A.ap(inverse(C), 3)
      #if sum(map(len, CC4)) < 50:
      #  sclCC4 = scl(CC4)
      #else:
      #  sclCC4 = 1e10
      sclCC = min(sclCC, sclCC2)
      traintrack_matrices_sep_domain(TROLLOPDIR, A, ell, v)
      lowerbound = trollop(None, ell, rank=rank, mat_comp=True, sep_domain=True, method='CIPT')
      if lowerbound > 0:
        nontrivial_lowerbound_count += 1
      print "Testing endo: ", A, " with scl=", sclCC, " and lower bound= ", lowerbound
      sys.stdout.flush()
      if abs(lowerbound - sclCC) < 0.01:
        good_examples_count += 1
        examples.append( (A, v, C, ell) )
        print (A, v, C, ell)
        sys.stdout.flush()
        break
  return (all_endos_count, \
          nonstupid_expanding_fix_hom_endos_count, \
          nontrivial_lowerbound_count, \
          good_examples_count, \
          examples)
  
  
  
  
  
  
  
def census_of_endomorphisms(word_len) :
  W = all_words_of_len_le(word_len, ['a','b'])
  LW = len(W)
  num_endos = 0
  num_stupid = 0
  num_with_fixed = 0
  num_with_fixed_and_hom_iden = 0
  hom_idens = []
  good_endos = []
  num_with_fixed_with_good_power = 0
  num_with_fixed_with_good_power_nonzero_rots = 0
  num_with_fixed_and_hom_iden_nonzero_rots = 0
  num_expanding = 0
  CO = cyclic_orders(2)
  for i in xrange(LW):
    for j in xrange(LW): 
      num_endos += 1
      if stupid_endomorphism((W[i], W[j])):
        num_stupid += 1
        continue 
      A = morph({'a':W[i], 'b':W[j]})
      H = homology_matrix(A)
      v = fixed_vector(A)      
      if v == None:
        continue
      if H == [[1,0],[0,1]]:
        num_with_fixed_and_hom_iden += 1
      num_with_fixed += 1
      if not A.is_expanding():
        continue
      num_expanding += 1
      AT = A.tripod_action()
      if AT==None:
        print "This shouldn't happen"
        continue
      good_power = 1
      good_orders = []
      while True:
        good_orders = [c for c in CO if AT.preserves_order(c, good_power)]
        if len(good_orders)>0 or good_power == 4:
          break
        else:
          good_power += 1
      if len(good_orders) == 0:
        continue 
      num_with_fixed_with_good_power += 1
      C = chain_from_vector(v)
      iC = cyc_red(A.ap(inverse(C), good_power))
      rots = [(c, rot(c, C+iC)) for c in good_orders]
      if any([x[-1]!=0 for x in rots]):
        num_with_fixed_with_good_power_nonzero_rots += 1
        if H == [[1,0],[0,1]]:
          num_with_fixed_and_hom_iden_nonzero_rots += 1
      
      A_entry = [A, H, v, rots]
      good_endos.append(A_entry)
      if H == [[1,0],[0,1]]:
        hom_idens.append(A_entry)
        
    print "\r", i, " of ", LW,
    sys.stdout.flush()
  print "\nFound:"
  print num_endos, " endos"
  print num_stupid, " stupid endos"
  print num_with_fixed, " with a fixed vector"
  print num_expanding, "ish of these are expanding"
  print num_with_fixed_with_good_power, " with a good power <= 4"
  print num_with_fixed_with_good_power_nonzero_rots, " with a good power and at least one nonzero rot"
  print num_with_fixed_and_hom_iden, " which fixed homology"
  print num_with_fixed_and_hom_iden_nonzero_rots, " of the homology fixers with a good rot"
  return (good_endos, hom_idens)


def find_good_geometric_endos(word_len, rank=2, trials=100000):
  if rank==2:
    W = all_words_of_len_le(word_len, alphabet[:rank])
    do_random = False
    LW = len(W)
    num_to_do = LW*LW
  else:
    do_random = True
    num_to_do = trials
    
  CO = cyclic_orders(rank)
  good_endos = []
  geometrics_found = 0
  chains_too_big = 0
  found = 0
  for i in xrange(num_to_do):
    print "\r",i," of ", num_to_do,
    print "; ", geometrics_found, " geometric; ",
    print chains_too_big, " chains too big; ",
    print found, " found",
    sys.stdout.flush()
    if do_random:
      targets = [ random_reduced_word(word_len, rank) for j in xrange(rank) ]
    else:
      targets = [ W[i/LW], W[i%LW] ]
    if stupid_endomorphism(targets):
      continue
    A = morph( dict( zip( alphabet[:rank], targets ) ) )
    v = fixed_vector(A)      
    if v == None:
      continue
    if not A.is_expanding():
      continue
    AT = A.tripod_action()
    if AT==None:
      print "This shouldn't happen"
      continue
    good_power = 1
    good_orders = []
    while True:
      #good_ordersm1 = [c for c in CO if AT.preserves_order(c, good_power-1)]
      good_orders = [c for c in CO if AT.preserves_order(c, good_power)]
      if len(good_orders)>0 or good_power == 4:
        break
      else:
        good_power += 1
    if len(good_orders) == 0:
      continue 
    geometrics_found += 1
    C = chain_from_vector(v)
    iC = cyc_red(A.ap(inverse(C), good_power))
    CC = C + iC
    #print "computing scl(", CC,")"; sys.stdout.flush()
    sCC = smart_scl(CC)
    if sCC == None:
      chains_too_big += 1
      continue
    #print "Done"; sys.stdout.flush()
    rots = [(c, rot(c, C+iC)) for c in good_orders]
    successful_orders = []
    for r in rots:
      if r[-1] == 2*sCC:
        successful_orders.append(r)
    if len(successful_orders) > 0:
      found += 1
      good_endos.append( (A, good_power, C, sCC, successful_orders) )
      print good_endos[-1]
  return good_endos



def find_good_counting_quasis(num_auts, \
                              num_words_per_aut, \
                              word_len, \
                              aut_let, \
                              chain_type='simple_word', \
                              quasi_len=2, \
                              num_scls=2, \
                              target_words=None, \
                              rank=2, \
                              scylla=False) :
  hom_list = []
  possible_surfaces = []
  try:
    ### Create the list of endomorphisms and fixed vectors
    if target_words==None:
      for i in xrange(num_auts):
        A, hom_image = random_homomorphism_with_fixed_vector(2, \
                                          random.randint(aut_len/2, aut_len), \
                                          random.randint(aut_len/2, aut_len))
        if cyc_red(A.ap('a')) == 'a' or cyc_red(A.ap('b')) == 'b':
          continue
        else:
          hom_list.append( (A,hom_image) )
    else:
      list_lens = [len(TW) for TW in target_words]
      all_target_indices = gen_tuples(list_lens)
      print "Created list of all target tuples"
      sys.stdout.flush()
      #print all_target_indices
      for TI in all_target_indices:
        w = [target_words[i][TI[i]] for i in xrange(len(target_words))]
        if stupid_endomorphism(w):
          continue
        A = morph(dict([(alphabet[i], w[i]) for i in xrange(len(w))]))
        ev = fixed_vector(A)
        if ev != None:
          hom_list.append( (A, ev) )
  
    print "Done creating hom list"
    sys.stdout.flush()
    #print hom_list
    ### now run through all the endomorphisms
    for i in xrange(len(hom_list)):
      chain_list = []
      A, hom_image = hom_list[i]
      possible_surfaces.append([A])
      print "Trying hom: ", A, " with fixed vector: ", hom_image
      
      ###make the test chain(s)
      if chain_type=='simple_word':
        chain_list = [ (hom_image[k]*alphabet[k] if hom_image[k] >=0 \
                       else (-hom_image[k])*inverse(alphabet[k]))      \
                            for k in xrange(len(hom_image)) ]
        chain_list = [[''.join(chain_list)]]
      elif chain_type=='simple_chain':
        chain_list = [ (hom_image[k]*alphabet[k] if hom_image[k] >=0 \
                       else (-hom_image[k])*inverse(alphabet[k]))      \
                     for k in xrange(len(hom_image)) ]
        chain_list = [chain_list]
      else:
        if random_word_with_hom_image(word_len, 2, hom_image) == None:
          chain_list = [random_word_with_hom_image(word_len-1, 2, hom_image) \
                                              for i in xrange(num_words_per_aut)]
        else:
          chain_list = [random_word_with_hom_image(word_len, 2, hom_image) \
                                              for i in xrange(num_words_per_aut)]
        chain_list = [chain_list]
      
      ### now run through the test chains
      for C in chain_list:
        print "Trying chain: ", C
        print "Length ", quasi_len, " quasi values: ",
        images = [A.iterate(i).ap(C) for i in xrange(1,5)]
        q_vals = [trollop(images[i] + inverse(C), quasi_len) for i in xrange(4)]
        print "quasi values, k=1..4: ", q_vals
        if q_vals[-1] > q_vals[0]:
          print "Computing scls: ", 
          scl_vals = [smart_scl(images[i] + inverse(C)) for i in xrange(num_scls)]
          print scl_vals
          possible_surfaces[-1].append([C, q_vals, scl_vals])
      
      if len(possible_surfaces[-1]) == 1:
        del possible_surfaces[-1]
      
    return possible_surfaces
  except None:
    return possible_surfaces
        
      











def find_surface_via_scl(num_auts, num_words_per_aut, aut_len, word_len, chain_type='words', target_words=None) :
  WH_gens = list(WAgen(2))
  possible_surfaces = []
  inc = 0
  count = 0
  aut_list = []
  if target_words == None:
    for i in xrange(num_auts):
      A, hom_image = random_homomorphism_with_fixed_vector(2, random.randint(aut_len/2, aut_len), random.randint(aut_len/2, aut_len))
      if cyc_red(A.ap('a')) == 'a' or cyc_red(A.ap('b')) == 'b':
        continue
      else:
        aut_list.append( (A,hom_image) )
  else:
    for i in xrange(1,len(target_words)):
      w1 = target_words[i]
      for j in xrange(1,len(target_words)):
        w2 = target_words[j]
        if w1 == w2 or \
           sorted([cyc_red(w1).lower(), cyc_red(w2).lower()]) == ['a','b'] or \
           cyc_red(w1).lower() == 'a' or cyc_red(w2).lower() == 'b':
          continue
        A = morph({'a':w1, 'b':w2})
        ev = fixed_vector(A)
        if ev != None:
          aut_list.append( (A, ev) )
  #print aut_list
  try:
    for i in xrange(len(aut_list)):
      A, hom_image = aut_list[i]
      print "Trying aut: ", A, " with fixed vector: ", hom_image
      
      for j in xrange(num_words_per_aut):
        if chain_type == 'words':
          w = [random_word_with_hom_image(word_len, 2, hom_image)]
          if w == [None]:
            w = [random_word_with_hom_image(word_len-1, 2, hom_image)]
        elif chain_type == 'basic_chains':
          w = [ (hom_image[k]*alphabet[k] if hom_image[k] >=0 \
                 else (-hom_image[k])*inverse(alphabet[k]))      \
                 for k in xrange(len(hom_image)) ]
          w = [''.join(w)]
          sys.stdout.flush()
        elif chain_type == 'chains':
          w = random_chain_with_hom_image(word_len, 2, word_len/3, hom_image)
          if w == None:
            w = random_chain_with_hom_image(word_len+1, 2, word_len/3, hom_image)
        else:
          w = []
        
        print "\tTrying chain: ", w, " scl1 = ",       
        sys.stdout.flush() 
        
        word_len = sum([len(x) for x in min_in_orbit(w + A.ap(inverse(w)))])
        print "Computing scl of ", min_in_orbit(w + A.ap(inverse(w)) )
        if word_len < 50:
          scl1 = scl( min_in_orbit(w + A.ap(inverse(w))) )
        elif word_len < 120:
          scl1 = matlab_scl(min_in_orbit(w + A.ap(inverse(w))))[0]
        else:
          print "Skipping"
          if chain_type == 'basic_chains':
            break
          continue
          
        print "got ", scl1, " scl2 = ",
        
        word_len = sum([len(x) for x in min_in_orbit(w + A.iterate(2).ap(inverse(w)))])
        if word_len < 50:
          scl2 = scl( min_in_orbit(w + A.iterate(2).ap(inverse(w))) )
        elif word_len < 120:
          scl2 = matlab_scl( min_in_orbit(w + A.iterate(2).ap(inverse(w))) )[0]
        else:
          print "Skipping"
          if chain_type == 'basic_chains':
            break
          continue
        
        if scl2 > 1.9*scl1:
          print scl2, ' ***********'
        else:
          print scl2
        sys.stdout.flush()
        count += 1
        inc += (scl2 - scl1)
        if (scl2 > 1.9*scl1) :
          print A, w, scl1, scl2
          word_len = sum([len(x) for x in min_in_orbit(w + A.iterate(3).ap(inverse(w)))])
          print "Checking next level -- word len = ", word_len
          sys.stdout.flush()
          if word_len > 150:
            print "Oops word len too long"
            possible_surfaces.append( [A, w, scl1, scl2] )
            if chain_type == 'basic_chains':
              break
            continue
          print "running matlab"
          try:
            scl3 = matlab_scl( min_in_orbit(w + A.iterate(3).ap(inverse(w))) )[0]
            print "Done -- scl3 = ", scl3
          except:
            print "matlab failed, continuing"
            scl3 = 0
          if scl3 > 0.9*3.0*scl1:
            print "passed!"
            possible_surfaces.append( [A, w, scl1, scl2, scl3] )
          else:
            print "failed"

        if chain_type=='basic_chains':
          break
          
    print "Average increase: ", float(inc)/float(count)
    return possible_surfaces
    
    
  except:
    print "Average increase: ", float(inc)/float(count)
    return possible_surfaces





