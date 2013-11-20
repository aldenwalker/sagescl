from sage.all import *
import word
import covering
import fatgraph
import ends



def find_extremal_transfer(C_in, max_degree=None, degree_list=None, verbose=1):
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
  
  if verbose > 1:
    print "Got scl = ", s
    if verbose > 2:
      print "Got extremal fatgraph: "
      print str(F)
  
  #these are the good ones
  found_transfers = []
  
  #go through the list of degrees
  for deg in deg_list:
    #go through all possible permuations for all the generators
    #i.e. all covering spaces (ouch)
    P = Permutations(deg).list()
    T = Tuples(P, rank)
    base_gens = word.alphabet[:rank]
    if verbose > 2:
      for t in T:
        G = covering.FISubgroup(base_gens, t)
        print "Found covering subgroup ", G
        GF = F.lift(G)
        print "Found lifted fatgraph: "
        print(GF)
        GFE = GF.ends()
        print "Found end list: ", GFE
        GFE_listed = [e.lift(G) for e in GFE]
        print "Found the lifted end list: ", GFE_lifted
        compat_orders = compatible_cyclic_orders(GFE_lifted, rank)
        print "Found compatible cyclic orders: ", compat_orders
        if len(compat_orders) != 0:
          print "*** good transfer ", (t, G, compat_orders)
          found_transfers.append((t, G, compat_orders))
    else:
      for t in T:
        G = covering.FISubgroup(base_gens, t)
        GF = F.lift(G)
        GFE = GF.ends()
        compat_orders = ends.compatible_cyclic_orders(GFE, rank)
        if len(compat_orders) != 0:
          if verbose > 1:
            print "*** good transfer ", (t, G, compat_orders)
          found_transfers.append((t, G, compat_orders))

  return found_transfers


















      
