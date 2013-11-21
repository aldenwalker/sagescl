import random as RAND

from sage.all import *

alphabet = list('abcdefghijklmnopqrstuvwxyz')

nextLetterChoices = {'a':['a','b','B'], 'b':['a','A','b'], 'A':['A','b','B'],'B':['a','A','B']}
L = ['a','b','c','A','B','C']
nextLetterChoices3 = dict( [ ( x, [y for y in L if x != y.swapcase()]) for x in L] )
L4 = ['a','b','c','d','A','B','C','D']
nextLetterChoices4 = dict( [ ( x, [y for y in L4 if x != y.swapcase()]) for x in L4] )

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
  """return a dictionary whose keys are letters and whose values are lists of indices"""
  ans = {}
  for i, let in enumerate(w):
    ans[let] = ans.get(let, []) + [i]
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
  if len(w) == 0 or len(w) == 1:
    return (w, '')
  i=0
  Lw = len(w)
  hLw = Lw/2
  while i < hLw and w[i] == w[-(i+1)].swapcase():
    i+=1
  return ( (w[i:-i], w[:i]) if i>0 else (w, '') )

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

  

    
def cyclic_subword(w, i, L):
  Lw = len(w)
  L_left = L
  ans = ''
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

def random_hom_triv_chain(n, rank=2):
  gens = alphabet[:rank]
  words = all_words_of_len(2, gens)
  all_gens = gens + inverse(gens)
  gen_indices = dict([(all_gens[i], i) for i in xrange(len(all_gens))])
  LW = len(words)
  ngens = len(all_gens)
  #choose the number of letters of each kind
  gen_counts = [0 for i in xrange(rank)]
  for i in xrange((n/2)+1):
    gen_counts[RAND.choice(xrange(rank))] += 1
  letter_vector = gen_counts + gen_counts #this counts all the letters
  #these record which positions are connected to where
  outgoing_positions = [ [None for j in xrange(letter_vector[i])] for i in xrange(ngens)]
  incoming_positions = [ [None for j in xrange(letter_vector[i])] for i in xrange(ngens)]
  #these record the indices that we can't connect
  inverse_positions = [all_gens.index(g.swapcase()) for g in all_gens]
  #these record the available positions -- this is slow
  available_outgoing_positions = [(i,j) for i in xrange(ngens) for j in xrange(letter_vector[i])]
  available_incoming_positions = [(i,j) for i in xrange(ngens) for j in xrange(letter_vector[i])]
  while len(available_outgoing_positions) > 0:
    #choose an available outgoing position
    k1 = RAND.choice(xrange(len(available_outgoing_positions)))
    (i1,j1) = available_outgoing_positions[k1]
    while True:
      lai = len(available_incoming_positions)
      k2 = RAND.choice(xrange(lai))
      (i2,j2) = available_incoming_positions[k2]
      if inverse_positions[i1] != i2:
        break
    outgoing_positions[i1][j1] = (i2, j2)
    incoming_positions[i2][j2] = (i1, j1)
    del available_outgoing_positions[k1]
    del available_incoming_positions[k2]
  #now we could choose a random permutation, but I think we don't need this?
  is_position_done = [ [False for j in xrange(letter_vector[i])] for i in xrange(ngens)]
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
    while pos != start_pos:
      output_word += all_gens[pos[0]]
      is_position_done[pos[0]][pos[1]] = True
      pos = outgoing_positions[pos[0]][pos[1]]
    chain.append(output_word)

  return chain
  
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
