from sage.all import *
import word
import covering
import fatgraph

#an end is a finite prefix, followed by an infinitely repeating word
class End:
  def __init__(self, prefix, repeat):
    self.prefix = prefix
    self.repeat = repeat
  
  def ap_morph(self, phi):
    ppre = phi.ap(self.prefix)
    prea = phi.ap(self.repeat)
    core, core_prefix = cyc_red_get_conjugate(prea)
    ppre = multiply_words(ppre, core_prefix)
    if ppre == '':
      return End('', core)
    while ppre[-1] == core[0].swapcase():
      ppre += core
    return End(ppre, core)
      



def find_extremal_transfer(C_in, max_degree=None, degree_list=None):
  """given a chain, try to find an extremal rot transfer.  This will 
  take a whole bunch of covers and then look through all *basic* rots 
  for all the covers"""
  
  #get a list of all the covering degrees we should go through
  if max_degee == None:
    if degree_list == None:
      print "I need a max degree or degree list"
      return None
    deg_list = degree_list_in
  else:
    deg_list = range(2,max_degree+1)
  
  #get the chain as a list of words
  if type(C_in) == str:
    C = [C_in]
  else:
    C = C_in
  
  #get the rank
  rank = word.chain_rank(C)
  
  #get the cwd
  cur_dir = os.getcwd()
  
  #get an extremal surface
  s = scl.scl(C, 'local', ['-o', cur_dir + '/temp_extremal_surface.fg'])
  s = Rational(s.numerator, s.denominator) #this turns it into a sage thing
  F = fatgraph.read_file(cur_dir + '/temp_extremal_surface.fg')
  
  #get the list of ends
  base_E = f.ends()
  
  #these are the good ones
  found_transfers = []
  
  #go through the list of degrees
  for deg in deg_list:
    #go through all possible permuations for all the generators
    #i.e. all covering spaces (ouch)
    P = Permutations(deg).list()
    T = Tuples(P, rank)
    base_gens = word.alphabet[:rank]
    for t in T:
      G = covering.FISubgroup(base_gens, t)
      conj_morphs = G.conjugation_action()
      lifted_ends = [[e.ap_morph(phi) for e in base_E] for phi in conj_morphs]
      compat_orders = cyclic_orders_compatible_with_end_lists(lifted_ends)
      if len(compat_orders) != 0:
        found_transfers.append((t, G, compat_orders))

  return found_transfers


















      