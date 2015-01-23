import random as RAND
import morph
import itertools

from sage.all import *

alphabet = list('abcdefghijklmnopqrstuvwxyz')

nextLetterChoices = {'a':['a','b','B'], 'b':['a','A','b'], 'A':['A','b','B'],'B':['a','A','B']}
L = ['a','b','c','A','B','C']
nextLetterChoices3 = dict( [ ( x, [y for y in L if x != y.swapcase()]) for x in L] )
L4 = ['a','b','c','d','A','B','C','D']
nextLetterChoices4 = dict( [ ( x, [y for y in L4 if x != y.swapcase()]) for x in L4] )

def next_letter_dict(rank):
  our_alphabet = alphabet[:rank] + [ell.swapcase() for ell in alphabet[:rank]]
  next_letters = dict([ (a,[b for b in our_alphabet if b!=a.swapcase()]) for a in our_alphabet])
  next_letters[''] = our_alphabet
  return next_letters

def edge_relative_to_id(pos, path):
  """gives the edge obtained by going to pos and then following path.
  If the last thing is to backtrack, it doesn't cancel that"""
  lp = len(path)
  current_word = pos
  for i,ell in enumerate(path):
    if len(current_word) == 0 or current_word[-1] != ell.swapcase():
      current_word += ell
    else:
      if i == lp-1:
        current_word += ell
      else:
        current_word = current_word[:-1]
  return current_word


def trim_trailing_cancel(w):
  if len(w)<2 or w[-1] != w[-2].swapcase():
    return w
  else:
    return w[:-1]


def non_reduced_index( w ):
  """return the first index i for which i, i+1 is nonreduced"""
  LW = len(w)
  if LW == 0 or LW == 1:
    return None
  for i in xrange(LW-1):
    if w[i] == w[i+1].swapcase():
      return i
  if w[LW-1] == w[0].swapcase():
    return LW-1
  return None

def letter_index_dict(w):
  """Return a dictionary whose keys are letters and whose values are lists of indices.
  If w is a list of words, then the values are lists of pairs (word_index, letter_index)."""
  if type(w) == list:
    sub_dicts = [letter_index_dict(W) for W in w]
    all_keys = list(set(itertools.chain(*[d.keys() for d in sub_dicts])))
    big_dict = {}
    for k in all_keys:
      for i, d in enumerate(sub_dicts):
        big_dict.setdefault(k, []).extend([(i,j) for j in d.get(k,[])])
    return big_dict
  ans = {}
  for i, let in enumerate(w):
    ans.setdefault(let, []).append(i)
  return ans


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

def multiply_words_marked(w1, w2=None) :
  if type(w1)==str:
    w1 = clean_markings(w1)
    w2 = clean_markings(w2)
    i=0
    j=0
    w1L = len(w1)
    w2L = len(w2)
    while i<w1L and j<w2L:
      if w1[w1L-i-1] == '.':
        i += 1
      if w2[j] == '.':
        j += 1
      if w1[w1L-i-1] != w2[j].swapcase():
        break
      i += 1
      j += 1
    return clean_markings(w1[:w1L-i] + '.' + w2[j:])
  elif type(w1) == list:
    return clean_markings('.' + reduce(multiply_words_marked, w1))

def cyclic_alignment(w1, i, w2, j):
  L = 0
  LW1 = len(w1)
  LW2 = len(w2)
  M = max(LW1, LW2)
  while L<M and w1[(i+L)%LW1] == w2[(j+L)%LW2].swapcase():
    L += 1
  return L

def rectangle_alignment(w1, i, w2, j):
  """returns the length of the alignment between w1 and w2
  when going *forward* along w1 and *backward* along w2 from positions 
  i and j, respectively, alignment here meaning the letters are inverse"""
  L = 0
  LW1 = len(w1)
  LW2 = len(w2)
  M = max(LW1, LW2)
  while L<M and w1[(i+L)%LW1] == w2[(j-L)%LW2].swapcase():
    L += 1
  return L

def clean_markings(w):
  i=0
  while i < len(w) and w[-1] == '.':
    w = w[-1] + w[0:-1]
    i += 1
  if i == len(w):
    return ''
  i=0
  while w[i] == '.':
    i += 1
  if i==0:
    return w
  else:
    return w[i-1:]

def cyc_red(w, marked=False):
  if type(w) == list:
    return [cyc_red(x, marked) for x in w]
  LW = len(w)
  if marked:
    if len(w) == 0:
      return w
    if len(w) == 1:
      if w == '.':
        return ''
      else:
        return w
    if w[0] == '.':
      return '.' + cyc_red(w[1:], True)
    
  if len(w) == 0 or len(w) == 1:
    return w
  else:
    i = 0
    while w[i] == w[LW-i-1].swapcase():
      i+=1
    return w[i:LW-i]

def cyc_red_get_conjugate(w):
  """Cyclically reduce, plus return the conjugating word"""
  if len(w) == 0 or len(w) == 1:
    return (w, '')
  i=0
  Lw = len(w)
  hLw = Lw/2
  while i < hLw and w[i] == w[-(i+1)].swapcase():
    i+=1
  return ( (w[i:-i], w[:i]) if i>0 else (w, '') )

def word_reduce(w):
  w2 = str(w)
  i=0
  while i < len(w2)-1:
    if w2[i] == w2[i+1].swapcase():
      w2 = w2[:i] + w2[i+2:]
      if i>0:
        i-=1
    else:
      i += 1
  return w2

def sign(letter):
  if letter.isupper():
    return -1
  else:
    return 1

def runs(w):
  wl = w.lower()
  s = wl.count('ab') + wl.count('ba')
  if wl[0] != wl[-1]:
    s += 1
  return s

def count_prefixes(w, p):
  """count the number of prefixes of p in w"""
  i=0
  n=0
  LP = len(p)
  while w[i:i+LP] == p:
    i += LP
    n += 1
  return n

    
def cyclic_subword(w, i, L):
  Lw = len(w)
  L_left = L
  ans = ''
  i = i%Lw
  while L_left > Lw:
    ans += w[i:] + w[:i]
    L_left -= Lw
  if i+L_left > Lw:
    return ans + w[i:] + w[:(i+L_left)%Lw]
  else:
    return w[i:i+L_left]


#return the subword which is including i1 through, not including, i2
def cyclic_subword_between_indices(w, i1_in, i2_in):
  Lw = len(w)
  i1 = i1_in%Lw
  i2 = i2_in%Lw
  ans = ''
  if i2 < i1:
    return w[i1:] + w[:i2]
  else:
    return w[i1:i2]



def inverse(w, marked=False):
  if type(w) == str:
    if marked and len(w)>0 and w[0] == '.':
      return '.' + w[:0:-1].swapcase()
    else :
      return w[::-1].swapcase()
  else:
    return [inverse(x,marked) for x in w]

def commutator(x,y):
  return multiply_words([x,y,inverse(x), inverse(y)])

def case_notation_single(w):
  if '^' not in w:
    return w
  gen, val = w.split('^')
  val = int(val)
  if val > 0:
    return val*gen
  else:
    return (-val)*gen.swapcase()
  
  
def is_inverse_chain(C1_in, C2_in):
  C1 = map(cyc_red, C1_in)
  C2 = map(cyc_red, C2_in)
  for c in C2:
    c = inverse(c)
    if not c in C1:
      return False
  for c in C1:
    c = inverse(c)
    if not c in C2:
      return False   
  return True
  
  
def power_reduce(w):
  Lw = len(w)
  for potential_period in xrange(1, len(w)+1):
    if Lw%potential_period != 0:
      continue
    p = Lw/potential_period
    if w == p*w[:potential_period]:
      return [w[:potential_period], p]
  #this will always return in the for loop
  return None

def min_cyclic(w):
  return min([w[i:] + w[:i] for i in xrange(len(w))])
  
def min_cyclic_get_shift(w):
  """returns the min cyclic word, plus the position in the original 
  word where the min cyclic one starts"""
  return min([ ( w[i:] + w[:i], i) for i in xrange(len(w))])

def cyclically_contained(x, L):
  return min_cyclic(x) in map(min_cyclic, L)
  
def last(x):
  return x[-1] 
  
def chain_rank(C):
  lets = ''.join(C)
  if len(lets) == 0:
    return 0
  ml = max(lets.lower())
  return ord(ml) - 97 + 1
 
def chain_rank_and_gens(C):
  """figures out exactly how many letters there are (not like the above function
  and what they are"""
  found_letters = set()
  lets = ''.join(C)
  for L in lets:
    found_letters.add(L)
  gens = set([x for x in found_letters])
  gens.union([x.swapcase() for x in found_letters])
  gens = [g for g in gens if g.islower()]
  return (len(gens), gens)
  
  
 
 
 
 
def pair_jump(a,b,w):
  posa,posb,posB = -1,-1,-1
  wl = len(w)
  B = inverse(b)
  if a == b or a == inverse(b):
    return 0
  for i in xrange(wl):
    if w[i] == a:
      posa=i
    elif w[i]==b:
      posb = i
    elif w[i]==B:
      posB = i
    if posa != -1 and posb != -1 and posB != -1:
      break
  if posa > posb:
    posb += wl
  if posa > posB:
    posB += wl
  if posb < posB:
    return posb-posa
  else:
    return -(wl-(posb-posa))
    
  

def rot(w, C_in):
  if type(C_in) == list:
    return sum([rot(w, x) for x in C_in])
  C = cyc_red(C_in[::-1])
  lc = len(C)
  ans = sum([ pair_jump(C[i], C[(i+1)%lc], w) for i in xrange(lc)])
  return ans/len(w)
  
  
def p2_rot(w_in, p, gen_2='b'):
  if type(w_in) == list:
    return sum([p2_rot(w,p,gen_2) for w in w_in])
  gen_1 = ('a' if gen_2=='b' else 'b')
  w = w_in
  if w[0] != gen_2:
    while w[0] != gen_2:
      w = w[-1] + w[:-1]
  groups = w.split(gen_2)[1:]
  return sum([ (1 if g==gen_1 else -1) for g in groups])
  
  
def all_words_of_len_le(n, gens):
  allGens = gens + [x.swapcase() for x in gens if x.swapcase() not in gens]
  words = [''] + allGens
  oldLen = 1
  for iters in xrange(n-1):
    startIndex = oldLen
    oldLen = len(words)
    for i in xrange(startIndex, oldLen):
      newWords = [words[i] + x for x in allGens if words[i][-1] != x.swapcase()]
      words.extend(newWords)
  return words

def all_words_of_len(n, gens):
  if n==0:
    return ['']
  allGens = gens + [x.swapcase() for x in gens if x.swapcase() not in gens]
  words = [''] + allGens
  oldLen = 1
  for iters in xrange(n-1):
    startIndex = oldLen
    oldLen = len(words)
    for i in xrange(startIndex, oldLen):
      newWords = [words[i] + x for x in allGens if words[i][-1] != x.swapcase()]
      words.extend(newWords)
  return words[oldLen:] 

def all_cyc_red_hom_triv_word_reps_of_len_le(n, gens):
  W = all_words_of_len_le(n, gens)
  S = set([least_cyclic_inverse_rep(w) for w in W if is_hom_triv(w) and w != ''])
  return S

def least_cyclic_rep(w):
  s = cyc_red(w)
  least_word = s
  for i in xrange(len(s)):
    s2 = s[i:] + s[:i]
    if s2 < least_word:
      least_word = s2
  return least_word
  
def least_cyclic_inverse_rep(w):
  return min( least_cyclic_rep(w), least_cyclic_rep(inverse(w)) )
  
def all_once_tagged_loops_of_len(n, gens):
  allGens = gens + [x.swapcase() for x in gens if x.swapcase() not in gens]
  words = [g for g in allGens]
  new_words = []
  for w in words:
    for g in allGens:
      if g != w[-1].swapcase():
        new_words.append(w + '/' + g + g.swapcase() + '/')
  words = new_words
  new_words = []
  for w in words:
    for g in allGens:
      if g != w[-2].swapcase() and g != w[0].swapcase():
        new_words.append(w + g)
  words = new_words
  
  oldLen = 0
  for iters in xrange(n-3):
    startIndex = oldLen
    oldLen = len(words)
    for i in xrange(startIndex, oldLen):
      newWords = [words[i] + x for x in allGens if words[i][-1] != x.swapcase()]
      words.extend(newWords)
  words = [w for w in words[oldLen:] if w[0] != w[-1].swapcase()]
  return words
  
def all_once_tagged_loop_extensions(f, L, gens):
  if len(f) == 0:
    return all_once_tagged_loops_of_len(L, gens)
  allGens = gens + [x.swapcase() for x in gens if x.swapcase() not in gens]
  words = [f + g for g in allGens if g != f[-1].swapcase()]
  new_words = []
  for w in words:
    for g in allGens:
      if g != w[-1].swapcase():
        new_words.append(w + '/' + g + g.swapcase() + '/')
  words = new_words
  new_words = []
  for w in words:
    for g in allGens:
      if g != w[-2].swapcase() and g != w[0].swapcase():
        new_words.append(w + g)
  words = new_words
  
  oldLen = 0
  while len(words[-1]) < L:
    startIndex = oldLen
    oldLen = len(words)
    for i in xrange(startIndex, oldLen):
      newWords = [words[i] + x for x in allGens if words[i][-1] != x.swapcase()]
      words.extend(newWords)
  words = [w for w in words[oldLen:] if w[0] != w[-1].swapcase()]
  return words
  
  
  
  
  
  
def vector_to_chain(v, word_list):
  c = []
  for i in xrange(len(v)):
    if v[i] != 0:
      c.append(str(v[i]) + word_list[i])
  return c
  
def random_reduced_word(n, rank=2):
  if rank == 2:
    letters = [RAND.choice(['a','A','b','B'])]
    for i in range(n-1):
      letters.append(RAND.choice(nextLetterChoices[letters[i]]))
    return ''.join(letters)
  elif rank == 3:
    letters = [RAND.choice(L)]
    for i in range(n-1):
      letters.append(RAND.choice(nextLetterChoices3[letters[i]]))
    return ''.join(letters)
  elif rank == 4:
    letters = [RAND.choice(L4)]
    for i in range(n-1):
      letters.append(RAND.choice(nextLetterChoices4[letters[i]]))
    return ''.join(letters)    
  else:
    all_letters = []
    next_choices = {}
    for i in xrange(rank):
      all_letters.extend([alphabet[i], inverse(alphabet[i])])
    for i in xrange(rank):
      next_choices[all_letters[2*i]] = all_letters[:2*i+1] + all_letters[2*i+2:]
      next_choices[all_letters[2*i+1]] = all_letters[:2*i] + all_letters[2*i+1:]
#print next_choices
    word = [RAND.choice(all_letters)]
    for i in xrange(n-1):
      word.append( RAND.choice(next_choices[word[-1]]) )
    return ''.join(word)

def simplify_finite_gen_power(power, order):
  if abs(power - order) < power:
    return power - order
  else:
    return power

def word_list_reduce(LC):
  L = [x for x in LC]
  for i in xrange(len(L)-1):
    if L[i] == '':
      continue
    left_word = i
    right_word = i+1
    while L[left_word][-1] == L[right_word][0].swapcase():
      L[left_word] = L[left_word][:-1]
      if L[left_word] == '':
        if left_word == 0:
          L[right_word] = L[right_word][1:]
          break
        else:
          left_word -= 1
      L[right_word] = L[right_word][1:]
      if L[right_word] == '':
        if right_word == len(L):
          break
        else:
          right_word += 1
  return L    

#returns the vector recording which (cyclic) subwords appear in w
def subword_vector(w, word_list):
  ww = cyc_red(w)
  ell = len(word_list[0])
  lww = len(ww)
  ww = ((ell/lww)+2)*ww
  vec = [0 for i in xrange(len(word_list))]
  for i in xrange(lww):
    word = ww[i:i+ell]
    vec[ word_list.index(word) ] += 1
  return vec

  
def parse_subword_vector(w, word_list):
  s = w.split(' ')
  v = [0 for word in word_list]
  for x in s:
    i = -1
    while -i <= len(x) and x[i].isalpha():
      i -= 1
    inp_word = x[i+1:]
    coef = x[:i+1]
    i = word_list.index(inp_word)
    v[i] += sage_eval(coef) # v[i] += int(coef)
  return v

def print_subword_vector(v, word_list):
  s = [str(v[i]) + word_list[i] for i in xrange(len(v)) if v[i] != 0]
  return ' '.join(s)

#find a match in the second coordinate in L1[start1:end1] and L2[start2:end2]
def find_match(L1, start1, end1, L2,start2, end2):
  for i in xrange(start1, end1):
    for j in xrange(start2, end2):
      if L1[i][1] == L2[j][1]:
        return (i,j)
  return None
  
#swap the elements at positions i and j
def list_swap(L, i, j): 
  temp = L[i]
  L[i] = L[j]
  L[j] = temp
 
#L1 and L2 are lists of elements of the form (x,n), where n is number
#this sorts them so that as many as possible of the n match up
#so [(a,1), (b,2), (c,2)], [(d,2), (e,3), (f,1)] -> puts the 2's and 1's together
def match_sort(L1, L2):
  sorted_index = 0
  LL1 = len(L1)
  LL2 = len(L2)
  while True:
    #find matching numbers
    match_values = find_match(L1, sorted_index, LL1, L2, sorted_index, LL2)
    if match_values == None:
      return
    i,j = match_values
    list_swap(L1, sorted_index, i)
    list_swap(L2, sorted_index, j)
    sorted_index += 1
    if sorted_index == min(LL1, LL2):
      return

#joins a path in a traintrack into a word
def coalesce(L):
  return ''.join([x[0] for x in L])

#produces a chain from a traintrack       
def subword_vector_to_chain(v, word_list):
  gluings = {}
  indices = dict( [ (word_list[i], i) for i in xrange(len(word_list))] )
  #examine each vertex
  verts = list(set([w[:-1] for w in word_list]))
  for vert in verts:
    #print "Working vertex ", vert
    incoming = [ [w, v[indices[w]]] for w in word_list if w[1:] == vert and v[indices[w]] != 0]
    outgoing = [ [w, v[indices[w]]] for w in word_list if w[:-1] == vert and v[indices[w]] != 0]
    match_sort(incoming, outgoing)
    #print "Got sorted edges: ", incoming, outgoing
    while len(incoming) > 0:
      if incoming[0][1] < outgoing[0][1]:
        gluings[(incoming[0][0], outgoing[0][0])] = incoming[0][1]
        outgoing[0][1] -= incoming[0][1]
        del(incoming[0])
      elif incoming[0][1] > outgoing[0][1]:
        gluings[(incoming[0][0], outgoing[0][0])] = outgoing[0][1]
        incoming[0][1] -= outgoing[0][1]
        del(outgoing[0])
      else:
        gluings[(incoming[0][0], outgoing[0][0])] = outgoing[0][1]
        del(outgoing[0])
        del(incoming[0])
      #print "Did one gluing: ", incoming, outgoing
  #print "Gluings: ", gluings
  #now go back and put the pieces together
  complete_paths = []
  remaining_v = [x for x in v]
  dim = len(remaining_v)
  while True:
    #print "Complete paths: ", complete_paths
    #print "Gluings left: ", gluings
    i=0
    while i < dim and remaining_v[i] == 0:
      i += 1
    if i == dim:
      break
    #build a loop out of this one
    current_path = [word_list[i]]
    coefficient = v[i]
    gluing_list = []
    #print "Starting with ", current_path
    while True:
      #print "Current path: ", current_path
      current_edge = current_path[-1]
      #do we close up the loop?
      for i in xrange(len(current_path)):
        close_up = gluings.get((current_edge, current_path[i]), 0)
        if close_up > 0:
          close_up_index = i
          break
      if close_up > 0:
        gluing_list.append( (current_edge, current_path[close_up_index]) )
        current_path = current_path[close_up_index:]
        gluing_list = gluing_list[close_up_index:]
        coefficient = min(coefficient, close_up)
        for g in gluing_list:
          gluings[g] -= coefficient
          remaining_v[indices[g[0]]] -= coefficient
          if gluings[g] == 0:
            del gluings[g]
        complete_paths.append( (coefficient, coalesce(current_path) ) )
        break
      else:
        #I don't think it matter which we choose
        all_potential_gluings = [g for g in gluings if g[0] == current_path[-1]]
        all_potential_gluings.sort(key=lambda x:-gluings[x])
        current_path.append(all_potential_gluings[0][1])
        coefficient = min(coefficient, gluings[all_potential_gluings[0]])
        gluing_list.append(all_potential_gluings[0])
  
  #sometimes, we get multiple copies of the same word, so let's collect those
  complete_paths = [(x[0], min_cyclic(x[1])) for x in complete_paths]
  i = 0
  while i < len(complete_paths)-1:
    j = i+1
    while j < len(complete_paths):
      if complete_paths[i][1] == complete_paths[j][1]:
        complete_paths[i] = (complete_paths[i][0] + complete_paths[j][0], complete_paths[i][1])
        del complete_paths[j]
      else:
        j += 1
    i += 1
    
  return complete_paths
        
  
def is_hom_triv(C_in, rank=None):
  C = ([C_in] if type(C_in)==str else C_in)
  r = (rank if rank != None else chain_rank(C))
  for let in alphabet[:r]:
    countl = sum([w.count(let) for w in C])
    countu = sum([w.count(let.swapcase()) for w in C])
    if countl != countu:
      return False
  return True
  

def is_diskbusting(w, rank=None):
  r = (rank if rank != None else chain_rank(w))
  WG = whitehead_graph(morph.min_in_orbit(w,rank), rank)
  return WG.is_connected()

def whitehead_graph(C_in, rank=None):
  C = ([C_in] if type(C_in) == str else C_in)
  rank = (rank if rank != None else chain_rank(C))
  WG = Graph()
  for i in xrange(rank):
    WG.add_vertex(alphabet[i])
    WG.add_vertex(alphabet[i].swapcase())
  for w in C:
    wL = len(w)
    for i in xrange(wL):
      WG.add_edge( ( w[i].swapcase(), w[(i+1)%wL] ) )
  return WG

    
def random_reduced_finite_word(n, orders, first=None):
  num_gens = len(orders)
  W = []
  for i in xrange(n):
    if i==0 and first!=None:
      gen = first
    else:
      gen = RAND.randint(0, num_gens-1)
    while len(W) > 0 and gen == W[-1][0]:
      gen = RAND.randint(0, num_gens-1)
    W.append( ( gen, simplify_finite_gen_power(RAND.randint(1, orders[gen]-1), orders[gen]) ))
  #print W
  W = [ (x[-1]*alphabet[x[0]] if x[-1] > 0 else (-x[-1])*inverse(alphabet[x[0]])) for x in W]
  return ''.join(W)

    

def random_hom_triv_word(n, rank=2):
  if n%2 != 0:
    return None
  if rank==2:
    w = random_reduced_word(n)
    while w.count('a') != w.count('A') or w.count('b') != w.count('B'):
      w = random_reduced_word(n)
    return w
  elif rank==3:
    w = random_reduced_word(n, rank)
    while w.count('a') != w.count('A') or w.count('b') != w.count('B') or w.count('c') != w.count('C'):
      w = random_reduced_word(n, rank)
    return w
  else:
    w = random_reduced_word(n, rank)
    while w.count('a') != w.count('A') or \
          w.count('b') != w.count('B') or \
          w.count('c') != w.count('C') or \
          w.count('d') != w.count('D'):
      w = random_reduced_word(n, rank)
    return w 
  
def random_hom_triv_cyc_red_word(n, rank=2):
  w = random_hom_triv_word(n,rank)
  while w[0] == w[-1].swapcase():
    w = random_hom_triv_word(n,rank)
  return w

def random_cyc_red_word(n, rank=2):
  w = random_reduced_word(n, rank)
  while w[0] == w[-1].swapcase():
    w = random_reduced_word(n, rank)
  return w


def random_thrice_punctured_sphere(n, rank=2):
  w1 = random_reduced_word(n, rank)
  w2 = random_reduced_word(n, rank)
  w1w2 = inverse(multiply_words(w1, w2))
  return [cyc_red(x) for x in [w1, w2, w1w2]]
  
  
def old_random_hom_triv_chain(n, maxWords, rank=2):
  wordLens = []
  s = 0
  for i in xrange(maxWords-1):
    wordLens.append(RAND.randint(1,int(n)))
    s += wordLens[-1]
    if s > n:
      wordLens[-1] -= s-n
      break
    elif s==n:
      break
  if s<n:
    wordLens.append(n-s)
  print wordLens
  words = [random_cyc_red_word(x, rank) for x in wordLens]
  i=0
  while not is_hom_triv(words) and i<5:
    words = [random_cyc_red_word(x, rank) for x in wordLens]
    i+=1
  if not is_hom_triv(words):
    while not is_hom_triv(words):
      words = [cyc_red(random_reduced_word(x, rank)) for x in wordLens]
  
  return words

def random_hom_triv_chain(n, Rank=2, Gens=None):
  
  if Gens == None:
    rank = Rank
    gen_translation = None
  else:
    rank = len(Gens)
    gen_translation = dict( [(alphabet[i], Gens[i]) for i in xrange(rank)] \
                           +[(alphabet[i].swapcase(), Gens[i].swapcase()) for i in xrange(rank)] ) 
  gens = alphabet[:rank]
    
  all_gens = gens + inverse(gens)
  gen_indices = dict([(all_gens[i], i) for i in xrange(len(all_gens))])
  ngens = len(all_gens)
  #get a random word
  W = random_hom_triv_cyc_red_word(n, rank)
  #reassemble it
  outgoing_positions = [ [] for i in xrange(ngens)]
  start_pos = (gen_indices[W[0]], 0)
  for i in xrange(len(W)-1):
    #make a spot for our current position
    outgoing_positions[gen_indices[W[i]]].append(None)
    target_gen_len = len(outgoing_positions[gen_indices[W[i+1]]])
    outgoing_positions[gen_indices[W[i]]][-1] = (gen_indices[W[i+1]], target_gen_len)  
  #add the last one
  outgoing_positions[gen_indices[W[-1]]].append(start_pos)
  #apply some perumations
  gen_perms = [ Permutations(range(len(outgoing_positions[i]))).random_element() for i in xrange(ngens)]
  #follow it around
  is_position_done = [ [False for j in xrange(len(outgoing_positions[i]))] for i in xrange(ngens)]
  chain = []
  while True:
    start_pos = None
    for i in xrange(ngens):
      try:
        j = is_position_done[i].index(False)
        start_pos = (i,j)
        break
      except:
        pass
    if start_pos == None:
      break
    output_word = all_gens[start_pos[0]]
    is_position_done[start_pos[0]][start_pos[1]] = True
    pos = outgoing_positions[start_pos[0]][start_pos[1]]
    pos = (pos[0], gen_perms[pos[0]][pos[1]])
    while pos != start_pos:
      #print "Current pos: ", pos
      #print "Start pos: ", start_pos
      #print "Current output word: ", output_word
      output_word += all_gens[pos[0]]
      is_position_done[pos[0]][pos[1]] = True
      pos = outgoing_positions[pos[0]][pos[1]]
      pos = (pos[0], gen_perms[pos[0]][pos[1]])
    chain.append(output_word)
    
  #print "Chain before messing: ", chain
  chain = scl_reduce(chain)
  if gen_translation != None:
    chain = [ ''.join([gen_translation[x] for x in c]) for c in chain]
  #print "Chain after messing: ", chain
  return chain

def scl_reduce(C):
  chain = [min_cyclic(cyc_red(x)) for x in C]
  chain = [power_reduce(x) for x in chain]
  chain_words = [x[0] for x in chain]
  chain_powers = [x[1] for x in chain]
  LC = len(chain)
  for i in xrange(LC):
    for j in xrange(i+1, LC):
      if chain_words[i] == chain_words[j]:
        chain_powers[i] += chain_powers[j]
        chain_powers[j] = 0
  #now there's no duplicates
  for i in xrange(LC):
    try:
      j = chain_words.index(inverse(chain_words[i]))
      if chain_powers[i] >= chain_powers[j]:
        chain_powers[i] -= chain_powers[j]
        chain_powers[j] = 0
    except:
      pass
  chain = [chain_powers[i]*chain_words[i] for i in xrange(len(chain_powers)) if chain_powers[i]>0]
  return chain


class FreeGroupFamily:
  def __init__(self, C, w1, i1, w2, i2):
    #print "Making family from ", C, w1, i1, w2, i2
    self.C = C
    self.w1 = w1
    self.i1 = i1
    self.w2 = w2
    self.i2 = i2
    
  def __repr__(self):
    if self.w1 != self.w2:
      W1 = self.C[self.w1][:self.i1] + self.C[self.w1][self.i1] + '^n' + self.C[self.w1][self.i1+1:]
      W2 = self.C[self.w2][:self.i2] + self.C[self.w2][self.i2] + '^n' + self.C[self.w2][self.i2+1:]
      return ' + '.join([self.C[i] for i in xrange(len(self.C)) if i != self.w1 and i != self.w2] + [W1] + [W2])
    else:
      i = min(self.i1, self.i2)
      I = max(self.i1, self.i2)
      W = self.C[self.w1][:i] + self.C[self.w1][i] + '^n' + self.C[self.w1][i+1:I] + self.C[self.w1][I] + '^n' + self.C[self.w1][I+1:]
      return ' + '.join([self.C[j] for j in xrange(len(self.C)) if j != self.w1] + [W])

  def __call__(self, n):
    if self.w1 != self.w2:
      W1 = self.C[self.w1][:self.i1] + n*self.C[self.w1][self.i1] + self.C[self.w1][self.i1+1:]
      W2 = self.C[self.w2][:self.i2] + n*self.C[self.w2][self.i2] + self.C[self.w2][self.i2+1:]
      return [self.C[i] for i in xrange(len(self.C)) if i != self.w1 and i != self.w2] + [W1] + [W2]
    else:
      i = min(self.i1, self.i2)
      I = max(self.i1, self.i2)
      W = self.C[self.w1][:i] + n*self.C[self.w1][i] + self.C[self.w1][i+1:I] + n*self.C[self.w1][I] + self.C[self.w1][I+1:]
      return [self.C[j] for j in xrange(len(self.C)) if j != self.w1] + [W]  


def random_family(n,rank=2, gens=None):
  C = random_hom_triv_chain(n,rank,gens)
  while len(C) == 0:
    C = random_hom_triv_chain(n,rank,gens)
  #print "Making random family from ", C
  LC = len(C)
  w1 = RAND.choice(xrange(LC))
  i1 = RAND.choice(xrange(len(C[w1])))
  w2 = RAND.choice(xrange(LC))
  i2 = RAND.choice(xrange(len(C[w2])))
  while C[w1][i1] != C[w2][i2].swapcase():
    w2 = RAND.choice(xrange(LC))
    i2 = RAND.choice(xrange(len(C[w2])))
  return FreeGroupFamily(C, w1, i1, w2, i2)
    



  
def random_word_with_hom_image(n, rank, hom_image) :
  if (n - sum([abs(x) for x in hom_image]) ) % 2 != 0:
    #print "That hom image is impossible with that length"
    return None
  gens = alphabet[:rank]
  while True:
    w = random_reduced_word(n, rank)
    if hom_image == [ w.count(x) - w.count(inverse(x)) for x in gens]:
      return w
  
def random_chain_with_hom_image(n, rank, max_words, hom_image) :
  if (n - sum([abs(x) for x in hom_image]) ) % 2 != 0:
    #print "That hom image is impossible with that length"
    return None
  wordLens = []
  s = 0
  gens = alphabet[:rank]
  for i in xrange(max_words-1):
    wordLens.append(random.randint(1,int(n/2)))
    s += wordLens[-1]
    if s > n:
      wordLens[-1] -= s-n
      break
    elif s==n:
      break
  if s<n:
    wordLens.append(n-s)
  #print wordLens
  while True:
    words = [cyc_red(random_reduced_word(x)) for x in wordLens]
    if hom_image == [ sum([w.count(x) - w.count(inverse(x)) for w in words]) for x in gens]:
      return words
   


def any_equal_pairs(t, inds, signs):
  lets = [(t[inds[i]] if signs[i] > 0 else inverse(t[inds[i]])) for i in xrange(len(inds))]
  if lets[0] == lets[1] or lets[0] == lets[2] or lets[1] == lets[2]:
    return True
  return False

def all_cubes(initial_edges=None):
  if initial_edges == None:
    T = Tuples(['a','b','A','B'],12)
  else:
    T = Tuples(['a','b','A','B'],12-len(initial_edges))
  good_T = []
  for i in xrange(T.cardinality()):
    t = initial_edges + T[i]
    if i%1000 == 0:
      print '\r' + str(i),
      sys.stdout.flush()
    #check all the vertices -- each label goes from 
    #the lower vertex to the higher one
    #vert 0: 0, 3, 4
    if any_equal_pairs(t, (0, 3, 4), (1,1,1) ):
      continue
    #vert 1: -0, 1, 5
    if any_equal_pairs(t, (0,1,5), (-1,1,1) ):
      continue
    #vert 2: -1, -2, 6
    if any_equal_pairs(t, (1,2,6), (-1,-1,1) ):
      continue
    #vert 3: 2, -3, 7
    if any_equal_pairs(t, (2,3,7), (1,-1,1) ):
      continue
    #vert 4: -4, 8, 11
    if any_equal_pairs(t, (4,8,11), (-1,1,1) ):
      continue
    #vert 5: -5,9, -8
    if any_equal_pairs(t, (5,9,8), (-1,1,-1) ):
      continue
    #vert 6: -9, -6, -10
    if any_equal_pairs(t, (9,6,10), (-1,-1,-1) ):
      continue
    #vert 7: -7, 10, -11
    if any_equal_pairs(t, (7,10,11), (-1,1,-1) ):
      continue
    good_T.append(t)
  return good_T



def TPS_with_given_tripods(T1, T2):
  #0->2, 1->1, 2->0
  #the middle words are outgoing from T1
  gens = list(set( (''.join(T1+T2)).lower() ))
  middle_words = []
  for i in xrange(3):
    i1 = i
    i2 = 2-i
    if T1[i1][-1].lower() == T2[i2][-1].lower(): 
      # just insert another generator
      middle_words.append( (gens[0] if gens[0] != T1[i1][-1].lower() else gens[1]) )
    else: 
      #just insert T1[i1][-1] again
      middle_words.append(T1[i1][-1])
  ans = []
  for i in xrange(3):
    i1 = i
    i2 = 2-i
    w = T1[i1] + middle_words[i1] + inverse(T2[i2]) + T2[(i2+1)%3] + inverse(middle_words[(i1-1)%3]) + inverse(T1[(i1-1)%3])
    ans.append(w)
  return ans

def tripods_are_negative(T1, T2):
  try:
    i = T2.index(T1[0])
  except ValueError:
    return False
  if T2[(i+1)%3] == T1[2] and T2[(i+2)%3] == T1[1]:
    return True
  return False

def tripods_are_equal(T1, T2):
  try:
    i = T2.index(T1[0])
  except ValueError:
    return False
  if T2[(i+1)%3] == T1[1] and T2[(i+2)%3] == T1[2]:
    return True
  return False


def find_all_rectangles(L_in):
  if type(L_in) == str:
    L = [L_in]
  else:
    L = L_in
  rectangles = []
  num_words = len(L)
  for i in xrange(num_words):
    LW1 = len(L[i])
    for j in xrange(num_words):
      LW2 = len(L[j])
      for posi in xrange(LW1):
        for posj in xrange(LW2):
          ca = rectangle_alignment(L[i], posi, L[j], posj-1)
          if ca > 0:
            rectangles.append( ((i,posi), (j, posj), ca) )
  return rectangles
  

def find_maximal_rectangles(L_in):
  if type(L_in)==str:
    L = [L_in]
  else:
    L = L_in
  num_words = len(L)
  rectangles = set()
  for i in xrange(num_words):
    Li = len(L[i])
    for j in xrange(i, num_words):
      Lj = len(L[j])
      for posi in xrange(Li):
        for posj in xrange(Lj):
          #print "Trying positions ", posi, " and ", posj
          forward_alignment = rectangle_alignment(L[i], posi, L[j], posj-1)
          backward_alignment = rectangle_alignment(L[j], posj, L[i], posi-1)
          #print "Got forward and backward alignments of ", forward_alignment, " and ", backward_alignment
          starti = (posi-backward_alignment)%Li
          startj = (posj+backward_alignment)%Lj
          #print "Adding: ", ( (i, starti), (j,startj), forward_alignment + backward_alignment )
          rectangles.add( ( (i, starti), (j,startj), forward_alignment + backward_alignment ) )
  return rectangles


def find_maximal_tripod(L_in, verbose=0):
  if type(L_in)==str:
    L = [L_in]
  else:
    L = L_in
  num_words = len(L)
  #for all triples of positions, find the rectangle alignments
  positions = [(i,j) for i in xrange(num_words) for j in xrange(len(L[i]))]
  num_pos = len(positions)
  triples = Tuples(positions, 3)
  max_len = 0
  best_tripod = None
  for [(i0,j0),(i1,j1),(i2,j2)] in triples:
    ra0 = rectangle_alignment(L[i0], j0, L[i1], j1-1)
    ra1 = rectangle_alignment(L[i1], j1, L[i2], j2-1)
    ra2 = rectangle_alignment(L[i2], j2, L[i0], j0-1)
    total_length = ra0 + ra1 + ra2
    if total_length > max_len:
      max_len = total_length
      best_tripod = [(i0,j0),(i1,j1),(i2,j2), (ra0, ra1, ra2)]
      if verbose>0:
        print "New max length of ", total_length, " with ", [(i0,j0),(i1,j1),(i2,j2)]
  return (max_len, best_tripod)

def find_maximal_quadpod(L_in, verbose=0):
  if type(L_in)==str:
    L = [L_in]
  else:
    L = L_in
  #find all rectangles
  R = find_all_rectangles(L)
  R_dict = dict([ ((I1, I2), ell) for (I1, I2, ell) in R])
  MR = find_maximal_rectangles(L)
  num_words = len(L)
  positions = [(i,j) for i in xrange(num_words) for j in xrange(len(L[i]))]
  max_length = 0
  best_quadpod = None
  for ( inside0_left, inside1_left, middle_length) in MR:
    inside0_right = (inside0_left[0], (inside0_left[1]+middle_length)%len(L[inside0_left[0]]) )
    inside1_right = (inside1_left[0], (inside1_left[1]-middle_length)%len(L[inside1_left[0]]) )
    #choose the outside positions for the left and right junctions
    for outside_left in positions:
      ra0_left = R_dict.get( (inside1_left, outside_left), 0 )
      ra1_left = R_dict.get( (outside_left, inside0_left), 0 )
      for outside_right in positions:
        ra0_right = R_dict.get( (inside0_right, outside_right), 0 )
        ra1_right = R_dict.get( (outside_right, inside1_right), 0 )
        total_length = middle_length + (0.5*(ra0_left + ra1_left + ra0_right + ra1_right))
        if total_length > max_length:
          max_length = total_length
          best_quadpod = (total_length, (inside0_left, inside1_left, outside_left, inside0_right, outside_right, inside1_right))
  return best_quadpod
        


def find_maximal_quadpod_old(L_in, verbose=0):
  if type(L_in)==str:
    L = [L_in]
  else:
    L = L_in
  num_words = len(L)
  #for all triples of positions, find the rectangle alignments
  positions = [(i,j) for i in xrange(num_words) for j in xrange(len(L[i]))]
  num_pos = len(positions)
  quadruples = Tuples(positions, 4)
  max_len = 0
  best_quadpod = None
  for [(i0,j0),(i1,j1),(i2,j2),(other_i1, other_j1)] in quadruples:
    #if L[i0][j0] == L[i2][j2-1].swapcase() or \
    #   L[i1][j1] == L[i0][j0-1].swapcase() or \
    #   L[i2][j2] == L[i1][j1-1].swapcase():
    if L[i1][j1] == L[i0][j0-1].swapcase():
      continue #it's not reduced
    ra0 = rectangle_alignment(L[i0], j0, L[i1], j1-1)
    ra1 = rectangle_alignment(L[i1], j1, L[i2], j2-1)
    ra2 = rectangle_alignment(L[i2], j2, L[i0], j0-1)
    #go to the other end of rectangle 0
    other_i0, other_j0, other_i2, other_j2 = i0, (j0+ra0)%len(L[i0]), i1, (j1-ra0)%len(L[i1])
    #if L[i0][j0] == L[i2][j2-1].swapcase() or \
    #   L[i1][j1] == L[i0][j0-1].swapcase() or \
    #   L[i2][j2] == L[i1][j1-1].swapcase():
    if L[other_i0][other_j0] == L[other_i2][other_j2-1].swapcase():
      continue #it's not reduced
    other_ra0 = rectangle_alignment(L[other_i0], other_j0, L[other_i1], other_j1-1)
    other_ra1 = rectangle_alignment(L[other_i1], other_j1, L[other_i2], other_j2-1)
    #this should be the same as rectangle 0
    if verbose>0:
      other_ra2 = rectangle_alignment(L[other_i2], other_j2, L[other_i0], other_j0-1)
      if other_ra2 != ra0:
        print "Rectangles don't agree?"
        print (i0,j0),(i1,j1),(i2,j2),(other_i0, other_j0), (other_i1, other_j1), (other_i2, other_j2)
    inside_length = ra0
    outside_length = ra1 + ra2 + other_ra0 + other_ra1   
    total_length = inside_length + 0.5*outside_length
    if total_length > max_len:
      max_len = total_length
      best_quadpod = [(i0,j0),(i1,j1),(i2,j2), (other_i0, other_j0), (other_i1, other_j1), (other_i2, other_j2), (inside_length, outside_length, total_length), (ra1, ra2, ra0, other_ra0, other_ra1)]
      if verbose>0:
        print "New max length of ", total_length, " with ", [(i0,j0),(i1,j1),(i2,j2)]
  return (max_len, best_quadpod)

def find_best_vertex_pair(L_in, verbose=0):
  if type(L_in)==str:
    L = [L_in]
  else:
    L = L_in
  #find all rectangles
  R = find_all_rectangles(L)
  R_dict = dict([ ((I1, I2), ell) for (I1, I2, ell) in R])
  MR = find_maximal_rectangles(L)
  num_words = len(L)
  positions = [(i,j) for i in xrange(num_words) for j in xrange(len(L[i]))]
  max_length = 0
  best_quadpod = None
  for ( inside0_left, inside1_left, middle_length) in MR:
    inside0_right = (inside0_left[0], (inside0_left[1]+middle_length)%len(L[inside0_left[0]]) )
    inside1_right = (inside1_left[0], (inside1_left[1]-middle_length)%len(L[inside1_left[0]]) )
    #choose the outside positions for the left and right junctions
    for outside_left in positions:
      ra0_left = R_dict.get( (inside1_left, outside_left), 0 )
      ra1_left = R_dict.get( (outside_left, inside0_left), 0 )
      for outside_right in positions:
        ra0_right = R_dict.get( (inside0_right, outside_right), 0 )
        ra1_right = R_dict.get( (outside_right, inside1_right), 0 )
        total_length = middle_length + (0.5*(ra0_left + ra1_left + ra0_right + ra1_right))
        if total_length > max_length:
          max_length = total_length
          best_quadpod = (total_length, (inside0_left, inside1_left, outside_left, inside0_right, outside_right, inside1_right))
  return best_quadpod


def subwords(w_in, ell):
  if type(w_in) == list:
    return list(itertools.chain(*[subwords(w,ell) for w in w_in]))
  return [cyclic_subword(w_in, i, ell) for i in xrange(len(w_in))]





def minimal_unique_word_length(C, w_ind, ind):
  """find the minimal length L such that the word C[w_ind] starting 
  at ind has the property that its inverse doesn't appear in C"""
  L = 1
  Lw = len(C[w_ind])
  current_ind = ind
  inverse_positions = [(i,j) for (i,j) in [(i,j) for i in xrange(len(C)) for j in xrange(len(C[i]))] if C[i][j] == C[w_ind][ind].swapcase()]
  while len(inverse_positions) > 0:
    #print inverse_positions
    current_ind = (current_ind+1)%Lw
    inverse_positions = [(i,(j-1)%len(C[i])) for (i,j) in inverse_positions]
    inverse_positions = [(i,j) for (i,j) in inverse_positions if C[i][j] == C[w_ind][current_ind].swapcase()]
    L += 1
  return L

def distinct_words(C_in, try_all_starts=False):
  """try to find the maximum number of words in C such that their inverse don't appear.
  It returns the number of them found.  If this number is x, note that it gives 
  a lower bound scl(C) >= x/12"""
  if type(C_in) == str:
    C = [C_in]
  else:
    C = C_in
  nw = len(C)
  found_subwords = []
  for i in xrange(nw):
    Lw = len(C[i])
    start_indices = ([0] if not try_all_starts else range(Lw))
    max_found_subwords = []
    for j in start_indices:
      current_start_index = j
      distance_to_start_index = Lw
      found_subwords_this_start_index = []
      while True:
        mwl = minimal_unique_word_length(C, i, current_start_index)
        distance_to_start_index -= mwl
        current_start_index = (current_start_index + mwl)%Lw
        if distance_to_start_index < 0:
          break
        else:
          found_subwords_this_start_index.append( (i,current_start_index, mwl) )
          if distance_to_start_index == 0:
            break
      if len(found_subwords_this_start_index) > len(max_found_subwords):
        max_found_subwords = found_subwords_this_start_index
    found_subwords += max_found_subwords
    
  return found_subwords


def reduce_using_relators(w, R, cyclically=False):
  """Reduce a word using the relators given.  This assumes a Dehn presentation
  (it searchs for locations where more than half a relator appears).  If 
  cyclically=True, it reduces it cyclically."""
  #make the mapping on words from the relators
  if type(w) == list:
    return [reduce_using_relators(W) for W in w]
  if type(R) == str:
    R = [R]
  R_sym = R + inverse(R)
  LR_sym = map(len, R_sym)
  bigger_than_half = [(ell/2)+1 for ell in LR_sym]
  letters_in_relators = letter_index_dict(R_sym)
  rw = w
  #print "R_sym: ", R_sym
  #print "letters_in_relators: ", letters_in_relators
  while True:
    #try to find a non-reduced spot
    Lrw = len(rw)
    non_reduced_spot = None
    for i in xrange(Lrw):
      agreements = letters_in_relators[rw[i]]
      agreement_len = 1
      good_agreements = []
      while len(agreements) > 0:
        if cyclically==False and i + agreement_len == Lrw:
          #check if there are any good agreements
          for a in agreements:
            if agreement_len >= bigger_than_half[a[0]]:
              good_agreements.append( (a[0], a[1], i, agreement_len) )
          break
        new_agreements = []
        #print "Agreements: ", agreements, " length: ", agreement_len
        for a in agreements:
          if LR_sym[a[0]] == agreement_len or \
             R_sym[a[0]][(a[1] + agreement_len)%LR_sym[a[0]]] != rw[(i+agreement_len)%Lrw]:
            #this agreement doesn't continue; if it's long enough, it's good, though
            if agreement_len >= bigger_than_half[a[0]]:
              good_agreements.append( (a[0], a[1], i, agreement_len) )
          else:
            new_agreements.append(a)
        agreements = new_agreements
        agreement_len += 1
      #print "Good agreeements: ", good_agreements
      if len(good_agreements) > 0:
        non_reduced_spot = max(good_agreements, key=lambda x:x[-1])
        break
    
    if non_reduced_spot == None:
      break
    #splice in the inverse of the relator
    ri, rij, rwi, al = non_reduced_spot
    r_remainder = cyclic_subword( R_sym[ri], rij+al, LR_sym[ri]-al )
    to_splice = inverse(r_remainder)
    #print "To splice: ", to_splice
    if cyclically == True:
      #rotate so we don't worry about it rolling off the end
      #our replacement now starts right at the beginning of the word
      rw = rw[rwi:] + rw[:rwi]
      rw = to_splice + rw[al:]
      rw = word_reduce(rw)
      rw = cyc_red(rw)
    else:
      #if it's not cyclically, we don't worry about it running over anyway
      rw = rw[:rwi] + to_splice + rw[(rwi+al):]
      rw = word_reduce(rw)
  
  return rw
      
        























  




