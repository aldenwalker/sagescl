from sage.all import *

from word import *
from scl import *

def random_L_word(n,p):
  w = random_reduced_finite_word(n,[p,2],first=1)
  return 'bA' + w + 'bA'

def find_stutters(n, p, trials, max_len=30, word_in=None, pr=False):
  st = []
  for i in xrange(trials):
    w = (random_L_word(n,p) if word_in==None else word_in)
    vals = []
    j=0
    while (2*j+len(w))<max_len  and (len(vals) < 6 or any([x[1]!=x[2] for x in vals[-4:]])):
      w_temp = (j*'ba') + w
      vals.append( (j, p2_rot(w_temp,p,'b')/Integer(12), scylla('a3b2',w_temp,method='GUROBI')) )
      j+=1
    if pr:
      print w, vals
    st_pattern = [('1' if x[1]==x[2] else '0') for x in vals]
    if '0' in st_pattern:
      st.append( (w, ''.join(st_pattern)) )
  
  #only look at stutters
  real_st = []
  for (w,s) in st:
    switches = [i for i in xrange(len(s)-1) if s[i]!=s[i+1] ]
    if len(switches) > 1:
      real_st.append((switches[-1]-switches[-2],w,s))
  
  return real_st

