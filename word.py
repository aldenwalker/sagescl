import random as RAND

from sage.all import *

alphabet = list('abcdefghijklmnopqrstuvwxyz')

nextLetterChoices = {'a':['a','b','B'], 'b':['a','A','b'], 'A':['A','b','B'],'B':['a','A','B']}
L = ['a','b','c','A','B','C']
nextLetterChoices3 = dict( [ ( x, [y for y in L if x != y.swapcase()]) for x in L] )
L4 = ['a','b','c','d','A','B','C','D']
nextLetterChoices4 = dict( [ ( x, [y for y in L4 if x != y.swapcase()]) for x in L4] )


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
    
def cyc_red(w):
  if type(w) == list:
    return [cyc_red(x) for x in w]
  LW = len(w)
  if len(w) == 0 or len(w) == 1:
    return w
  else:
    i = 0
    while w[i] == w[LW-i-1].swapcase():
      i+=1
    return w[i:LW-i]

def sign(letter):
  if letter.isupper():
    return -1
  else:
    return 1
    


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
  for potential_period in xrange(1, len(w)):
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
        
  
def is_hom_triv(C_in):
  C = ([C_in] if type(C_in)==str else C_in)
  r = chain_rank(C)
  for let in alphabet[:r]:
    countl = sum([w.count(let) for w in C])
    countu = sum([w.count(let.swapcase()) for w in C])
    if countl != countu:
      return False
  return True
  
    
def random_reduced_finite_word(n, orders):
  num_gens = len(orders)
  W = []
  for i in xrange(n):
    gen = random.randint(0, num_gens-1)
    while len(W) > 0 and gen == W[-1][0]:
      gen = random.randint(0, num_gens-1)
    W.append( ( gen, simplify_finite_gen_power(random.randint(1, orders[gen]-1), orders[gen]) ))
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

  
  
def random_hom_triv_chain(n, maxWords, rank=2):
  wordLens = []
  s = 0
  for i in xrange(maxWords-1):
    wordLens.append(2*random.randint(1,int(n/2)))
    s += wordLens[-1]
    if s > n:
      wordLens[-1] -= s-n
      break
    elif s==n:
      break
  if s<n:
    wordLens.append(n-s)
  #print wordLens
  words = [cyc_red(random_reduced_word(x)) for x in wordLens]
  while sum([w.count('a') for w in words]) != sum([w.count('A') for w in words]) \
        or sum([w.count('b') for w in words]) != sum([w.count('B') for w in words]):
    words = [cyc_red(random_reduced_word(x)) for x in wordLens]
  
  return words
  
  
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
   


  
