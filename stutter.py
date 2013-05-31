from sage.all import *

from word import *
from scl import *

def random_L_word(n,p):
  w = random_reduced_finite_word(n,[p,2],first=1)
  #rotate so a maximal string of ba's is first
  while w[1]=='A':
    w = w[2:] + w[:2]
  while w[-1] == 'a':
    w = w[-2:] + w[:-2]
  #cut off all but one of the first string
  while w[3] == 'a':
    w = w[2:]
  return w

#returns an L word which begins with a maximal cyclic run of "ba"'s
#e.g. bababAbabA
#such that rot is extremal
def random_L_word_rot_extremal(n,p):
  while True:
    w = random_reduced_finite_word(n,[p,2],first=1)
    r = p2_rot(w, p, 'b')/Integer(12)
    s = scylla('a3b2',w,method='GUROBI')
    if r == s:
      #rotate so a maximal string of ba's is first
      while w[1]=='A':
        w = w[2:] + w[:2]
      while w[-1] == 'a':
        w = w[-2:] + w[:-2]
      return (w,r)


def free_stutter_pattern(w, rot_order, word_to_prepend, max_pow, pr=False):
  vals = []
  j=0
  sc_gen_str = ''.join([ alphabet[i]+'0' for i in xrange(chain_rank(w))] )
  while j<max_pow and (len(vals) < 6 or any([x[1]!=x[2] for x in vals[-4:]])):
    w_temp = (j*word_to_prepend) + w
    vals.append( (j,rot(rot_order, w_temp)/Integer(2), scylla(sc_gen_str, w_temp, method='GUROBI')) )
    j+=1
  if pr:
    print w, vals
  st_pattern = [('1' if x[1]==x[2] else '0') for x in vals]
  return ''.join(st_pattern)

def stutter_pattern(w, p, max_len, pr=False):
  vals = []
  j=0
  while j<max_len  and (len(vals) < 6 or any([x[1]!=x[2] for x in vals[-4:]])):
    w_temp = (j*'ba') + w
    vals.append( (j, p2_rot(w_temp,p,'b')/Integer(12), scylla('a3b2',w_temp,method='GUROBI')) )
    j+=1
  if pr:
    print w, vals
  st_pattern = [('1' if x[1]==x[2] else '0') for x in vals]
  return ''.join(st_pattern)

def stutter_length(st):
  switches = [i for i in xrange(len(st)-1) if st[i]!=st[i+1] ]
  if len(switches) > 1:
    return switches[-1]-switches[-2]
  else:
    return 0

def find_stutters(n, p, trials, max_len=30, word_in=None, pr=False):
  st = []
  for i in xrange(trials):
    w = (random_L_word(n,p) if word_in==None else word_in)
    st_pattern = stutter_pattern(w, p, (max_len-len(w))/2, pr=pr )
    if '0' in st_pattern:
      st.append( (w, st_pattern) )
  
  #only look at stutters
  real_st = []
  for (w,s) in st:
    sl = stutter_length(s)
    if sl > 0:
      real_st.append((sl,w,s))
  
  return real_st


