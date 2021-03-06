#!/usr/bin/python

import fractions
import subprocess
import random as RAND
import sys
import os
import re

from morph import *
#load('morph.py')
from scl import *
from word import *
from rot_Sn_action import *
from sage.all import *

import covering
import fatgraph

alphabet = list('abcdefghijklmnopqrstuvwxyz')

def cont_frac(x, nterms):
  L = []
  X = x
  for i in xrange(nterms):
    L.append(floor(X))
    d = math.modf(X)[0]
    if abs(d) < 0.0000000000001:
      return L
    X = 1/d
  return L



def approx_rat(x, tol=0.00000000001):
  L = cont_frac(x, 20)
  L_frac = convergents(L)
  for i in xrange(len(L_frac)):
    if abs(L_frac[i]-x) < tol:
      return L_frac[i]
  return None


def lcm(a,b):
  return (a*b)/fractions.gcd(a,b)
  
def integerize_fractions_sage(L):
  denoms = [x.denominator() for x in L]
  g = reduce(lcm, denoms)
  cleared = [int(g*x) for x in L]
  g = reduce(fractions.gcd, cleared)
  return [x/g for x in cleared]
  
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
  


def cyclic_order_boundary(o):
  return boundary(o)

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
  
  if any([power_reduce(targets[k])[0].lower() == alphabet[k] for k in xrange(len(targets))]):
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
  letter_sign = lambda x: (1 if x.islower() else -1)
  hom_matrix = [ [(letter_sign(x[1]) if x[1].lower() == g else 0) for x in words] for g in alphabet[:rank] ]
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

#puts it in gap form, except puts a p after every letter
def to_gap_prime(s):
  gf = to_gap(s)
  return re.sub(r"([a-z])", r"\1p", gf)

def from_gap(s):
  ss = s.split('*')
  return ''.join(case_notation_single(w) for w in ss)   
  
def from_gap_prime(s):
  return from_gap(s.replace('p',''))

#remove as many t's as possible by killing off twT=A.ap(w)
def simplify_stable(A, w):
  t_ind = w.find('t')
  if t_ind == -1:
    return w
  T_ind = w.find('T',t_ind)
  if T_ind == -1:
    return w
  closer_t = w.find('t', t_ind+1, T_ind)
  while closer_t != -1:
    t_ind = closer_t
    closer_t = w.find('t', t_ind+1, T_ind)
  w_new = multiply_words([ w[:t_ind], A.ap(w[t_ind+1:T_ind]), w[T_ind+1:] ])
  return simplify_stable(A, w_new) 
  
#conjugate and stuff to remove instances of t
#here we take twT = A.ap(w)
def remove_stable_letter(A, stab, w):
  t_count = w.count('t') - w.count('T')
  t_count_stab = stab.count('t') - stab.count('T')
  if t_count_stab < 0:
    return 1/0
  if t_count > 0:
    new_w = multiply_words([w] + (t_count/t_count_stab)*[inverse(stab)])
  elif t_count < 0:
    new_w = multiply_words((-t_count/t_count_stab)*[stab] + [w])
  else:
    new_w = w 
  new_w = simplify_stable(A, new_w)
  #now we have zero signed t count
  #but we must conjugate if necessary to remove Txt
  t_count = new_w.count('t')
  m = ceil(Integer(t_count)/Integer(t_count_stab))
  new_w = multiply_words(m*[stab] + [new_w] + m*[inverse(stab)])
  return simplify_stable(A, new_w)
    


#use gap to find finite index subgroups
def gap_finite_index_subgroups_with_homology(A, index_bound, iteratively=False, find_HNN_structure=False):
  rank = len(A.rules)/2
  gens = alphabet[:rank]
  relators = [to_gap(multiply_words([inverse(g), 'T', A.rules[g], 't'])) for g in gens]
  gap.eval('LoadPackage("fga");')
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
  if not find_HNN_structure:
    return subgroups
  
  #find the HNN structure on each group.
  gap.eval('FS := FreeGroup(' + ','.join(['"' + g + 'p"' for g in gens]) + ');')
  gap.eval('AssignGeneratorVariables(FS);')
  subgroups_with_structures = []
  for S in subgroups:
    #first, we need to eliminate the t's in all generators except one 
    #which will be the stable letter
    ind = S[0]
    S_gens = gap.eval('GeneratorsOfGroup(subgroups[' + str(ind+1) + ']);')
    S_gens = [s.replace('[','').replace(']','').strip() for s in S_gens.split(',')]
    S_gens = [from_gap(s) for s in S_gens]
    t_counts = [s.count('t')-s.count('T') for s in S_gens]
    num_gens = len(S_gens)
    for i in xrange(num_gens):
      if t_counts[i] < 0:
        t_counts[i] *= -1
        S_gens[i] = inverse(S_gens[i])
    g = gcd(t_counts)
    if not g in t_counts:
      print "gcd isn't in there -- error"
      return 1/0
    stable_ind = t_counts.index(g)
    S_gens[stable_ind] = simplify_stable(A, S_gens[stable_ind])
    ###print ind
    ###print S_gens
    ###print t_counts
    ###print stable_ind
    for i in xrange(num_gens):
      if i == stable_ind:
        continue
      S_gens[i] = remove_stable_letter(A, S_gens[stable_ind], S_gens[i])
    ###print S_gens
    
    #now we have all the generators
    #we need to keep htting them with the stable letters until
    #they stay in the subgroup
    trivial = False
    current_gens = [S_gens[i] for i in xrange(num_gens) if i != stable_ind]
    stab = S_gens[stable_ind]
    gap_gens = [to_gap_prime(s) for s in current_gens]
    if '' in gap_gens:
      trivial = True
      subgroups_with_structures.append(list(S) + ['trivial'])
      continue
    gap.eval('S := Subgroup(FS, [' + ','.join(gap_gens) + ']);')
    while True:
      acted_on_gens = [simplify_stable(A, multiply_words([stab, s, inverse(stab)])) for s in current_gens]
      gap_acted_on_gens = [to_gap_prime(s) for s in acted_on_gens]
      if '' in gap_acted_on_gens:
        trivial = True
        break
      any_missing = False
      for i in xrange(len(acted_on_gens)):
        if not eval(gap.eval(gap_acted_on_gens[i] + ' in S;')):
          current_gens.append(acted_on_gens[i])
          gap_gens.append(gap_acted_on_gens[i])
          any_missing = True
      if not any_missing:
        break
      gap.eval('S := Subgroup(FS, [' + ','.join(gap_gens) + ']);')
    if trivial:
      subgroups_with_structures.append(list(S) + ['trivial'])
      continue
    ###print current_gens
    free_S_gens = gap.eval('mgs := MinimalGeneratingSet(S);')
    ###print free_S_gens
    free_S_gens = [s.translate(None, '[] \n') for s in free_S_gens.split(',')]
    free_S_gens = [from_gap_prime(s) for s in free_S_gens]
    ###print free_S_gens
    
    #now we have a list of the complete free generators
    #we need to figure out the endomorphism

    free_S_gens_targets = [simplify_stable(A, multiply_words([stab, s, inverse(stab)])) for s in free_S_gens]
    gap_free_S_gens_targets = [to_gap_prime(s) for s in free_S_gens_targets]
    free_S_rank = len(free_S_gens)
    gap.eval('H := FreeGroup( ' + str(free_S_rank) + ');')
    gap.eval('StoH := GroupHomomorphismByImages(S, H, mgs, GeneratorsOfGroup(H));')
    gap.eval('HtoS := GroupHomomorphismByImages(H, S, GeneratorsOfGroup(H), mgs);')
    gap.eval('Sphi := GroupHomomorphismByImages(S,S,mgs, [' +','.join(gap_free_S_gens_targets) + ']);')
    ###print gap.eval('HtoS*Sphi*StoH;')
    free_S_gens_targets = gap.eval('HtoS*Sphi*StoH;').split('->')[-1]
    #replace all the f's with letters
    while True:
      mo = re.search("f[0-9]+", free_S_gens_targets)
      if not mo:
        break
      val = int(free_S_gens_targets[mo.start()+1:mo.end()])
      free_S_gens_targets = free_S_gens_targets[:mo.start()] + \
                            alphabet[val-1] +                  \
                            free_S_gens_targets[mo.end():]
    free_S_gens_targets = [s.translate(None, ' \n][') for s in free_S_gens_targets.split(',')]
    free_S_gens_targets = [from_gap(s) for s in free_S_gens_targets]
    ###print free_S_gens_targets
    subgroup_A = morph(dict(zip(alphabet[:free_S_rank], free_S_gens_targets)))
    subgroups_with_structures.append( list(S[1:]) + [free_S_rank, subgroup_A] )

  return subgroups_with_structures
    



    
                          
  

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


def check_HNN_extensions_for_homology_covers(rank, word_len, index_check, endos=None, only_check_expanding=False, it=True, find_free_ranks=False):
  W = all_words_of_len(word_len, alphabet[:rank])
  LW = len(W)
  if endos == None:
    T = tuples_gen([LW for i in xrange(rank)])
  else:
    T = endos
  all_endos_count = 0
  nonstupid_expanding_count = 0
  trivially_good_count = 0
  good_count = 0
  examples = []
  non_examples = []
  for t in T:
    all_endos_count += 1
    if endos == None:
      targets = [W[i] for i in t]
      if stupid_endomorphism(targets):
        continue
      A = morph( dict( [(alphabet[i], targets[i]) for i in xrange(rank)]) )
    else:
      A = t
    if only_check_expanding and (not A.is_expanding()):
      continue
    nonstupid_expanding_count += 1
    print "Checking ", A, ' ...',
    sys.stdout.flush()
    if not find_free_ranks:
      subgroups_with_homology = gap_finite_index_subgroups_with_homology(A, index_check, iteratively=it)
    else:
      subgroups_with_homology = gap_finite_index_subgroups_with_homology(A, index_check, iteratively=it, find_HNN_structure=True)
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
        

def random_index_weighted_by_size(L):
  total_len = sum(map(len, L))
  r = RAND.random()
  s = len(L[0])*1.0/total_len
  i = 0
  while s < r:
    i += 1
    s += len(L[i])*1.0/total_len
  return i
    

#given an index, adjust it so that it's not on a tag
def move_off_tag(w, i, direction):
  Lw = len(w)
  if direction == 'forward':
    if w[i] == '/':
      if w[(i+1)%Lw] == w[(i+2)%Lw].swapcase():
        return (i+4)%len(w)
      else:
        return (i+1)%Lw
    elif w[i] == w[(i+1)%Lw].swapcase():
      return (i+3)%Lw
    elif w[i] == w[(i-1)%Lw].swapcase():
      return (i+2)%len(w)
  else:
    if w[i] == '/':
      if w[(i+1)%Lw] == w[(i+2)%Lw].swapcase():
        return (i-1)%len(w)
      else:
        return (i-4)%Lw
    elif w[i] == w[(i+1)%Lw].swapcase():
      return (i-2)%Lw
    elif w[i] == w[(i-1)%Lw].swapcase():
      return (i-3)%len(w)
  return i

#extend a matching as far as possible
#if folded, it'll return a zero length matching if it results in 
#a nonfolded surface
#also, if the folding ends right on a tag, it'll return zero length
#rather than a multi-tagged surface
def extend_matching(w1, i1_in, w2, i2_in, folded=False):
  if w1[i1_in] != w2[i2_in].swapcase():
    return (i1_in, 0, i2_in, 0, 0)
  L1 = 1
  L2 = 1
  Lw1 = len(w1)
  Lw2 = len(w2)
  i1 = i1_in
  i2 = i2_in
  total_match_len = 1
  #extend forward
  while True:
    tag1 = None
    if w1[(i1+L1)%Lw1] == '/':
      tag1 = w1[(i1+L1+1)%Lw1]
      L1 += 4
    tag2 = None
    if w2[(i2-1)%Lw2] == '/':
      tag2 = w2[(i2-3)%Lw2]
      i2 = (i2-4)%Lw2
      L2 += 4
    if w1[(i1+L1)%Lw1] == '.':
      L1 += 1
    if w2[(i2-1)%Lw2] == '.':
      i2 = (i2-1)%Lw2
      L2 += 1
      
    if w1[(i1+L1)%Lw1] != w2[(i2-1)%Lw2].swapcase():
      break
    if folded and tag1 != None and tag2 != None and tag1 == tag2:
      break
    total_match_len += 1
    L1 += 1
    i2 = (i2-1)%Lw2
    L2 += 1
    if L1 > Lw1 or L2 > Lw2:
      return (i1_in, 0, i2_in, 0, 0)
  #if we stopped right on a tag, don't use it
  #or if we stopped right by a mark
  if w1[(i1+L1-1)%Lw1] == '/' or w2[i2] == '/':
    return (i1_in, 0, i2_in, 0, 0)
  if w1[(i1+L1-1)%Lw1] == '.' or w2[(i2)%Lw2] == '.':
    return (i1_in, 0, i2_in, 0, 0)
  
  #print "Extended backwards to ", i1, L1, i2, L2
  
  #extend backward
  while True:
    tag1 = None
    if w1[(i1-1)%Lw1] == '/':
      tag1 = w1[(i1-3)%Lw1]
      i1 = (i1-4)%Lw1
      L1 += 4
    tag2 = None
    if w2[(i2+L2)%Lw2] == '/':
      tag2 = w2[(i2+L2+1)%Lw2]
      L2 += 4
    if w1[(i1-1)%Lw1] =='.':
      i1 = (i1-1)%Lw1
      L1 += 1
    if w2[(i2+L2)%Lw2] == '.':
      L2 += 1
    if w1[(i1-1)%Lw1] != w2[(i2+L2)%Lw2].swapcase():
      break
    if folded and tag1 != None and tag2 != None and tag1 == tag2:
      break
    total_match_len += 1
    i1 = (i1-1)%Lw1
    L1 += 1
    L2 += 1
    if L1 > Lw1 or L2 > Lw2:
      return (i1_in, 0, i2_in, 0, 0)
  #if we stopped right on a tag, don't use it
  if w1[i1] == '/' or w2[(i2+L2-1)%Lw2] == '/':
    return (i1_in, 0, i2_in, 0, 0)
  if w1[(i1)%Lw1] == '.' or w2[(i2+L2-1)%Lw2] == '.':
    return (i1_in, 0, i2_in, 0, 0)
  
  if L1 >= Lw1-1 or L2 >= Lw2-1:
    return (i1_in, 0, i2_in, 0, 0)
  return (i1, L1, i2, L2, total_match_len) 


#produce the chain which is obtained by gluing the specified words
def glue_tagged_chain(C, w1, i1, L1, w2, i2, L2):
  if w1 == w2:
    L = len(C[w1])
    #it's the same word
    first_loop = cyclic_subword_between_indices(C[w1],(i1+L1)%L, i2)
    tw1 = '/' + C[w1][i2] + C[w1][(i1+L1-1)%L] + '/' + first_loop
    second_loop = cyclic_subword_between_indices(C[w1],(i2+L2)%L, i1)
    tw2 = '/' + C[w1][i1] + C[w1][(i2+L2-1)%L] + '/' + second_loop
    return [tw1, tw2]
  else:
    Lw1 = len(C[w1])
    Lw2 = len(C[w2])
    #it's two different words
    return \
    ['/' + C[w2][i2] + C[w1][(i1+L1-1)%Lw1] + '/' + \
    cyclic_subword_between_indices(C[w1], (i1+L1)%Lw1, i1) + \
    '/' + C[w1][i1] + C[w2][(i2+L2-1)%Lw2] + '/' + \
    cyclic_subword_between_indices(C[w2], (i2+L2)%Lw2, i2)]



#find a random place to glue which glues up a lot of the letters involved
#returns a tagged chain
def random_folding(C_in):
  C = sorted([w for w in C_in], key=len, reverse=True)
  if len(C[0]) <= 10:
    return C
  word1 = random_index_weighted_by_size(C)
  while len(C[word1]) < 10:
    word1 = random_index_weighted_by_size(C)
  word2 = random_index_weighted_by_size(C)
  while len(C[word2]) < 10:
    word2 = random_index_weighted_by_size(C)
  best_gluing_length = 0
  best_gluing_indices = None
  #print "Chose words ", word1, word2, C[word1], C[word2]
  for i in xrange(tagged_len(C_in)):
    ind1 = RAND.randint(0, len(C[word1])-1)
    ind2 = RAND.randint(0, len(C[word2])-1)
    ind1 = move_off_tag(C[word1], ind1, 'forward')
    if C[word1][ind1] == '.':
      ind1 += 1
    ind2 = move_off_tag(C[word2], ind2, 'backward')   
    if C[word2][ind2] == '.':
      ind2 += 1
    #print "Trying matching", ind1, ind2
    (i1, l1, i2, l2, untagged_len) = extend_matching(C[word1], ind1, C[word2], ind2, folded=True)
    if untagged_len > 0 and \
       (C[word1][i1] != C[word2][(i2+l2-1)%len(C[word2])].swapcase() or \
       C[word1][(i1+l1-1)%len(C[word1])] != C[word2][i2].swapcase()):
      return 1/0
    #print "Found possible matching ", (i1, l1, i2, l2, untagged_len)
    if untagged_len > best_gluing_length:
      best_gluing_length = untagged_len
      best_gluing_indices = (i1, l1, i2, l2)
  # pair them
  if best_gluing_indices == None:
    return C
  ind1, len1, ind2, len2 = best_gluing_indices
  #print "Using matching ", (ind1, len1, ind2, len2)
  new_part = glue_tagged_chain(C, word1, ind1, len1, word2, ind2, len2)
  if any(['//' in x or 'A.a' in x or 'a.A' in x for x in new_part]):
    return 1/0
  return [C[i] for i in xrange(len(C)) if i not in (word1,word2)] + new_part
  



def check_all_HNN_extensions_for_ffolded(rank, word_len, power=1, endos=None, chains_to_try=None, time_limit=0, do_folding=False):
  TL = time_limit
  W = all_words_of_len(word_len, alphabet[:rank])
  LW = len(W)
  if endos == None:
    T = tuples_gen([LW for i in xrange(rank)]) 
  else:
    T = endos
  ns_endo_count = 0
  ffolded_count = 0
  found_one = False
  good_list = []
  bad_list = []
  torus_list = []
  for t in T:
    if endos==None:
      targets = [W[i] for i in t]
      if stupid_endomorphism(targets):
        continue
      A = morph( dict( [(alphabet[i], targets[i]) for i in xrange(rank)]) )
      if not A.is_expanding():
        continue
    else:
      A = t
    ns_endo_count += 1
    if chains_to_try == None:
      CL = [[random_hom_triv_cyc_red_word(8)] for i in xrange(4)]
    else:
      CL = chains_to_try #[['abAB'],['BA','a','b']]
    print "Trying ", A
    found_one = False
    for C in CL:
      if is_inverse_chain(C, inverse(A.ap(C))):
        print "torus"
        found_one = True
        torus_list.append( (A,C) )
        break
      IM = inverse(A.iterate(power).ap(C,marked=True), marked=True)
      IM = cyc_red(IM, marked=True)
      #print "Chain image: ", IM
      size = tagged_len(IM)
      print "chain image size ", size
      while do_folding:
        old_size = size
        IM = random_folding(IM)
        size = tagged_len(IM)
        if old_size - size < 10:
          break
        #print IM
      #print "Folded: ", IM
      print "Reduced to size ", tagged_len(IM)
      if tagged_len(IM) > 400:
        continue
      CC = C + IM
      print "With chain: ", C, ", ", CC
      if gallop('rose' + str(rank) + '.fg', CC, only_check_exists=True, folded=True, ffolded=len(C), solver="gurobi", time_limit=TL):
        print "It's good!\n\n\n"
        ffolded_count += 1
        good_list.append( (A,C) )
        found_one = True
        break
    if not found_one:
      print "No good"
      bad_list.append( A )
  return (ns_endo_count, torus_list, good_list, bad_list)
  

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
    if gallop('rose' + str(rank) + '.fg', C, only_check_exists=True, folded=True, ffolded=1, solver='gurobi'):
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
    
    

#this does not require expanding
#but it does iterate over orders
def find_extremal_surfaces(rank, word_len, endos=None, \
                           check_flat=True, check_counting=True, \
                           counting_ell=3):

  all_endos_count=0
  nonstupid_expanding_fix_hom_count = 0
  extremal_count=0
  flat_count=0
  counting_count=0
  good_endos=[]
  bad_endos=[]
  if endos==None:
    W = all_words_of_len(word_len, alphabet[:rank])
    LW = len(W)
    T=tuples_gen([LW for i in xrange(rank)])
    E = []
    for t in T:
      targets = [W[i] for i in t]
      all_endos_count += 1
      if stupid_endomorphism(targets):
        continue
      A = morph( dict( zip( alphabet[:rank], targets )))
      if not A.is_expanding() or len(A.fixed_space()) == 0:
        continue
      nonstupid_expanding_fix_hom_count += 1
      E.append(A)
  else:
    E = endos
    all_endos_count = len(E)
    E = [A for A in E if A.is_expanding() and len(A.fixed_space()) > 0]
    nonstupid_expanding_fix_hom_count = len(E)

  if check_flat:
    CO = cyclic_orders(rank)
    rots = [CQ(O) for O in CO]
    rots = [(Integer(1)/r.defect())*r for r in rots]

  for A in E:
    v = A.fixed_space()[0]
    C = [x for x in chain_from_vector(v) if x != '']
    iC = cyc_red(A.ap(inverse(C)))
    CC = C + iC
    print "For ", A, " with C=", C
    CCscl = gallop('rose'+str(rank)+'.fg', CC, solver="gurobi")
    CCscl = Rational(str(CCscl))
    print "scl=", CCscl
    found_extremal = False
    good_rots = []
    good_lower_bound = 0
    if check_flat:
      AT = tripod_action(A)
      preserved_order_inds = [i for i in xrange(len(CO)) if AT.preserves_order(CO[i])]
      rot_values = [(i, rots[i].ap(CC)) for i in preserved_order_inds]
      matched_rot_values = [x for x in rot_values if x[1] == 2*CCscl]
      print "Matched rot values: ", matched_rot_values
      if len(matched_rot_values) > 0:
        extremal_count += 1
        flat_count += 1
        is_extremal = True
        good_rots = matched_rot_values
        print "Good"
    
    if check_counting: # and not found_extremal:
      traintrack_matrices_sep_domain(TROLLOPDIR, A, counting_ell, v)
      lower_bound = trollop(None, counting_ell, rank=rank, mat_comp=True, sep_domain=True, method='GUROBI')
      print "Found counting lower bound: ", lower_bound
      lower_bound = Rational(str(lower_bound))
      #if lower_bound != CCscl:
      #  print "I find that ", lower_bound, " and ", CCscl, " arne't equal."
      if lower_bound == CCscl:
        if not found_extremal:
          extremal_count += 1
          found_extremal = True
        counting_count += 1
        good_lower_bound = lower_bound
        print "Good"
      
    if found_extremal:
      good_endos.append( (A, C, good_rots, good_lower_bound) )
    else:
      bad_endos.append( A )

  return (all_endos_count, nonstupid_expanding_fix_hom_count, extremal_count, \
          flat_count, counting_count, good_endos, bad_endos)

  
  

  
def find_flat_surfaces(rank, word_len, max_good_power=4):
  W = all_words_of_len(word_len, alphabet[:rank])
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
      if len(good_orders)>0 or good_power == max_good_power:
        break
      else:
        good_power += 1
    if len(good_orders) == 0:
      continue 
    good_power_count += 1
    if good_power > 2:
      continue
    C = [x for x in chain_from_vector(v) if x != '']
    iC = cyc_red(A.ap(inverse(C), good_power))
    rots = [(c, rot(c, C+iC)) for c in good_orders if rot(c, C+iC) > 0]
    if len(rots) == 0:
      continue
    nontrivial_rots_count += 1
    s = gallop('rose'+str(rank)+'.fg', C+iC, solver="gurobi")
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


def contains_all_gens(w, rank=2):
  if rank==2:
    return w.count('a') > 0 and w.count('A') > 0 and w.count('b') > 0 and w.count('B') > 0
  else:
    return False


def contains_aA(w):
  L = len(w)
  first_a = -1
  for i in xrange(L):
    if w[i] == 'a':
      first_a = i
      break
  if first_a == -1:
    return False
  last_A = -1
  for i in xrange(L-1, -1, -1):
    if w[i] == 'A':
      last_A = i
      break
  if last_A == -1 or last_A < first_a:
    return False
  return True

def fraction_that_satisfy(L, rank, test, trials=None):
  if trials == None:
    W = all_words_of_len(L, list(alphabet[:rank]))
    count = 0
    for w in W:
      if test(w):
        count += 1
    return count*1.0/len(W) 
  else:
    count = 0
    for i in xrange(trials):
      w = random_reduced_word(L, rank)
      if test(w):
        count += 1
    return count*1.0/trials
  

def matching_data(len1, len2):
  part1 = partitions_restricted(len1-4, range(len1-3), 5)
  part2 = partitions_restricted(len2-4, range(len2-3), 5)
  GS = {}
  for p1 in part1:
    arrange1 = arrangements(part1, 5)
    for p2 in part2:
      arrange2 = arrangements(part2, 5)
      for a1 in arrange1:
        for a2 in arrange2:
          compute_gluing_stats(a1, a2, GS)
          


def triples_match(t1, t2):
  if '.' in t1 or '.' in t2:
    return False
  return t1[1] == t2[1].swapcase() and t1[0] != t2[2].swapcase() and t1[2] != t2[0].swapcase()

def matches_all_triples(w, matching_data=None):
  triples = all_words_of_len(3, ['a','b'])
  matched = dict( zip( triples, len(triples)*[0]) )
  s = 0
  target = len(matched)
  if matching_data == None:
    matching_data = {}
    for t in triples:
      matching_data[t] = []
      for t2 in triples: 
        if triples_match(t, t2):
          matching_data[t].append(t2)
  
  for i in xrange(len(w)-2):
    t = w[i:i+3]
    for tm in matching_data[t]:
      if matched[tm] == 0:
        matched[tm] = 1
        s += 1
  
  return s == target


def triple_word_stats(L, trials):
  good_count = 0
  both_sides_good_count = 0
  bad_count = 0
  matching_data = None
  tL = 2*L+1
  for i in xrange(trials):
    w = random_reduced_word(2*L+1)
    if matches_all_triples(w, matching_data):
      good_count += 1
      if matches_all_triples(w[:L], matching_data) and \
         matches_all_triples(w[-L:], matching_data):
        both_sides_good_count += 1
    else:
      bad_count += 1
  
  print "Fraction good: ", good_count*1.0/trials
  print "Fraction really good: ", both_sides_good_count*1.0/trials
  print "Fraction bad: ", bad_count*1.0/trials
  print "After pairing, probability of segment:"
  print "Good but not really good (not paired): ", (good_count - both_sides_good_count)*1.0/trials
  print "Really good, not paired: ", (both_sides_good_count - bad_count)*1.0/trials
  print "Paired: ", 2*bad_count*1.0/trials


def random_good_not_really_good_word(L, matching_data):
  while True:
    w = random_reduced_word(2*L+1)
    if matches_all_triples(w, matching_data):
      if not matches_all_triples(w[:L], matching_data) or \
         not matches_all_triples(w[-L:], matching_data):
        return w

def random_really_good_word(L, matching_data):
  while True:
    w = random_reduced_word(2*L+1)
    if matches_all_triples(w[:L], matching_data) and \
       matches_all_triples(w[-L:], matching_data):
      return w

def random_bad_word(L, matching_data):
  while True:
    w = random_reduced_word(2*L+1)
    if not matches_all_triples(w, matching_data):
      return w

def random_paired(L, matching_data):
  first_is_good = ( RAND.random() < 0.5 )
  really_good_word = random_really_good_word(L, matching_data)
  bad_word = random_bad_word(L, matching_data)
  while not triples_match(really_good_word[L-1:L+2], bad_word[L-1:L+2]):
    bad_word = random_bad_word(L, matching_data)
  if first_is_good:
    return ( (really_good_word, bad_word), (True, False) )
  else :
    return ( (bad_word, really_good_word), (False, True) )

def make_tagged_word(w1, w2, loc1, loc2):
  return w1[:loc1] + '.' + w1[loc1] + w2[loc2] + '.' + w2[loc2+1:]

def select_word_from_probs(L, good, really_good, paired, matching_data):
  r = RAND.random()
  if r < good:
    return random_good_not_really_good_word(L, matching_data)
  elif r < good + really_good:
    return random_really_good_word(L, matching_data)
  else:
    p = random_paired(L, matching_data)
    return make_tagged_word(p[0][0], p[0][1], L, L)

def first_matching_triples(w1, w2):
  n = 0
  scan_len = min(len(w1), len(w2))-2
  while n < scan_len:
    for i in xrange(n):
      for j in xrange(n):
        if triples_match(w1[i:i+3], w2[-4-j:-1-j]):
          return (i,-4-j+len(w2))
    n += 1
  return None

#pair the beginning of tw1 with the *end* of tw2
def pinch_off_ends(tw1, tw2):
  #print "Pinching off ends from: ", tw1, tw2
  I = first_matching_triples(tw1, tw2)
  #print "Found index: ", I
  if I == None:
    return None
  word_ind_1, word_ind_2 = I
  word_ind_1 += 1
  word_ind_2 += 1
  t1 = make_tagged_word(tw1, tw2, word_ind_1, word_ind_2)
  #print "First word: ", t1
  W = make_tagged_word(tw2, tw1, word_ind_2, word_ind_1)
  t2, W = pinch(W)
  return t1, W, t2
  
  
#pinch off a tagged word (starting from the beginning)
#and return both words
def pinch(tw):
  L = floor(len(tw)/2)
  while L < len(tw) and (tw[L] == tw[L+1].swapcase() or tw[L] == tw[L-1].swapcase() or tw[L]=='.'):
    L+=1
  w1 = tw[:L]
  w2 = tw[L:]
  I = first_matching_triples(w1, w2)
  #print "Pinching off ", w1, w2
  if I == None:
    return [tw]
  else:
    i1 = I[0] + 1
    i2 = I[1] + 1
    return (w1[:i1] + '.' + w1[i1] + w2[i2] + '.' + w2[i2+1:], \
            w2[:i2] + '.' + w2[i2] + w1[i1] + '.' + w1[i1+1:])
  

def tagged_len(w):
  if type(w) == list:
    return sum(map(tagged_len, w))
  nt = w.count('/')/2
  return len(w) - 3*nt

def triple_gluing_stats(L, good, really_good, paired, trials):
  
  matching_data = None
  
  #these record the data about the words we create
  num_created = 0
  total_len = 0
  first_len = 0
  tag_len = 0
  
  circles_dict = {}
  tags_dict = {}
  
  middle_len = floor(L/2)-1
  for i in xrange(trials):
    #create the four words
    bottom = ['','']
    top = ['','']
    bottom[0] = select_word_from_probs(L, good, really_good, paired, matching_data)
    last_letter = bottom[0][-1]
    while True:
      bottom[1] = select_word_from_probs(L, good, really_good, paired, matching_data)
      first_letter = bottom[1][0]
      if first_letter != last_letter.swapcase():
        break
    top[0] = select_word_from_probs(L, good, really_good, paired, matching_data)
    last_letter = top[0][-1]
    while True:
      top[1] = select_word_from_probs(L, good, really_good, paired, matching_data)
      first_letter = top[1][0]
      if first_letter != last_letter.swapcase():
        break
    
    #print "Four words produced: "
    #print bottom, top
    
    bottom = bottom[0] + bottom[1] # join them
    top = top[0] + top[1]
    
    #print "Words: ", bottom, top
    
    #now pair the words recursively
    #pair the beginning and the end, and recurse
    t1, W, t2 = pinch_off_ends(bottom, top)
    #print "Pinched: ", t1, W, t2
    TL1 = tagged_len(t1)
    TL2 = tagged_len(t2)
    tag_len += TL1 + TL2
    tags_dict[TL1] = tags_dict.get(TL1, 0) + 1
    tags_dict[TL2] = tags_dict.get(TL2, 0) + 1
    first_len += tagged_len(W)
    #print "t1 has tagged len: ", tagged_len(t1)
    
    W2 = W
    while True:
      p = pinch(W2)
      #print "After pinching: ", p
      if len(p) == 1:
        num_created += 1
        TLp = tagged_len(p[0])
        circles_dict[TLp] = circles_dict.get(TLp, 0) + 1
        total_len += TLp
        break
      else:
        W1, W2 = p
        num_created += 1
        TLw = tagged_len(W1)
        total_len += TLw
        circles_dict[TLw] = circles_dict.get(TLw,0) + 1
  
  print "Total number of words created: ", num_created
  print "Average length: ", total_len*1.0/num_created
  print "Average first/last tag lengths: ", tag_len*1.0/(2*trials)  
  
  total_tags = sum(tags_dict.values())
  total_circles = sum(circles_dict.values())
  
  for t in tags_dict:
    tags_dict[t] = tags_dict[t]*1.0/total_tags
  for c in circles_dict:
    circles_dict[c] = circles_dict[c]*1.0/total_circles
  
  print "Fraction of words which are tags: "
  print "Tags length distribution: ", tags_dict
  print "Circles length distribution: ", circles_dict



#write out the transpose of V to the LP matrix, and put X on the RHS
def member_of_convex_cone(X, V, s="GLPK", interior=False):
  p = MixedIntegerLinearProgram(maximization=False, solver=s)
  #p.solver_parameter("simplex_or_intopt", "simplex_only")
  #p.solver_parameter("presolve_simplex", "GLP_ON")
  w = p.new_variable(real=True)
  for i in xrange(len(V[0])):
    new_eqn = 0
    nontrivial_eqn = False
    for j in xrange(len(V)):
      if V[j][i] != 0:
        nontrivial_eqn = True
        new_eqn += V[j][i]*w[j]
    if nontrivial_eqn:
      p.add_constraint(new_eqn == X[i])
  p.set_objective(0*w[0])
  if interior:
    for i in xrange(len(V)):
      p.set_min(w[i], 0.00001)
  try:
    soln = p.solve()
    return True
  except:
    return False


def find_nbhd_of_uniform_tagged(L):
  W = all_once_tagged_loops_of_len(L, ['a','b'])
  dim = len(W)
  hom_counts = [[w.count(g) - w.count(inverse(g)) for w in W] for g in ['a','b']]
  #go through and find a lot of small-density homologically trivial collections
  #small density means total length less than 40
  small_density_number = floor(60/L)
  if small_density_number*(L-1) % 2 != 0:
    small_density_number += 1
  print "Choosing small density collections of size", small_density_number
  #we expect to need something like len(W) of them to span
  good_vectors = []
  num_to_add = floor(1.2*dim)
  while True:
    for i in xrange(num_to_add):
      print "Trying to find new collection number", i
      while True: #try to find a good vector
        while True: #try to find a homologically trivial vector
          vec = [0 for j in xrange(dim)]
          vec_hom_counts = [0,0]
          for j in xrange(small_density_number):
            ind = RAND.randint(0, dim-1)
            vec[ind] += 1
            vec_hom_counts[0] += hom_counts[0][ind]
            vec_hom_counts[1] += hom_counts[1][ind]
          if vec_hom_counts == [0,0]:
            break
        #now we've got a homologically trivial vector; test if it bounds a folded surface
        word_vec = vector_to_chain(vec, W)
        if gallop('rose2.fg', word_vec, folded=True, only_check_exists=True):
          print "Found the good homologically trivial collection: ", ' '.join(word_vec)
          good_vectors.append(vec)
          break
    
    print "We have found ", len(good_vectors), " collections."
    print "Testing if the uniform vector is in the cone"
    if member_of_convex_cone([1 for i in xrange(dim)], good_vectors, s="Gurobi"):
      print "Yes; success"
      return good_vectors
      break
    print "No; failure"
    if len(good_vectors) > 1000000:
      print "Too many good vectors; bailing out"
      break
    num_to_add = floor(0.2*dim)
    print "Trying again; adding ", num_to_add, " collections"
  




def find_nbhd_of_uniform_tagged_fixed_portion(L, rank=2, num_to_fix=0):
  gens = list(alphabet[:rank])
  fixed_words = all_words_of_len(num_to_fix, gens)
  done_words = []
  good_collections = {}
  for f in fixed_words:
    if f in done_words:
      continue
    print "Doing fixed word", f
    done_words.extend([f,inverse(f)])
    W = all_once_tagged_loop_extensions(f, L, gens)
    dim = len(W)
    hom_counts = [[w.count(g) - w.count(inverse(g)) for w in W] for g in gens]
    #go through and find a lot of small-density homologically trivial collections
    #small density means total length less than something arbitrary
    small_density_number = floor(60/L)
    if small_density_number*(L-1) % 2 != 0:
      small_density_number += 1
    print "Choosing small density collections of size", small_density_number
    #we expect to need something like len(W) of them to span
    good_vectors = []
    num_to_add = floor(1.2*dim)
    success = False
    while True:
      for i in xrange(num_to_add):
        print "Trying to find new collection number", i
        while True: #try to find a good vector
          while True: #try to find a homologically trivial vector
            vec = [0 for j in xrange(dim)]
            vec_hom_counts = [0,0]
            for j in xrange(small_density_number):
              ind = RAND.randint(0, dim-1)
              vec[ind] += 1
              vec_hom_counts[0] += hom_counts[0][ind]
              vec_hom_counts[1] += hom_counts[1][ind]
            if vec_hom_counts == [0,0]:
              break
          #now we've got a homologically trivial vector; test if it bounds a folded surface
          word_vec = vector_to_chain(vec, W)
          print "Trying chain", word_vec
          if gallop('rose' + str(rank) + '.fg', word_vec, folded=True, only_check_exists=True):
            print "Found the good homologically trivial collection: ", ' '.join(word_vec)
            good_vectors.append(vec)
            break
      
      print "We have found ", len(good_vectors), " collections."
      print "Testing if the uniform vector is in the cone"
      if member_of_convex_cone([1 for i in xrange(dim)], good_vectors, s="Gurobi"):
        print "Yes; success"
        success = True
        break
      print "No; failure"
      if len(good_vectors) > 1000000:
        print "Too many good vectors; bailing out"
        break
      num_to_add = floor(0.2*dim)
      print "Trying again; adding ", num_to_add, " collections"
    if success:
      good_collections[f] = (W, good_vectors)
    else:
      print "We failed to find a collection for word", f
      return {}
    
  return good_collections
  
  
  





















