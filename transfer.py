from sage.all import *
import word
import covering
import fatgraph
import ends
import scl



def find_extremal_transfer(C_in, max_degree=None, degree_list=None, verbose=1):
  """given a chain, try to find an extremal rot transfer.  This will 
  take a whole bunch of covers and then look through all *basic* rots 
  for all the covers"""
  
  #get a list of all the covering degrees we should go through
  if max_degree == None:
    if degree_list == None:
      print "I need a max degree or degree list"
      return None
    deg_list = degree_list
  else:
    deg_list = range(2,max_degree+1)
  
  #get the chain as a list of words
  if type(C_in) == str:
    C = [C_in]
  else:
    C = C_in
  
  #get the rank
  rank, base_gens = word.chain_rank_and_gens(C)
  
  #get the cwd
  cur_dir = os.getcwd()
  
  #get an extremal surface
  s = scl.scl(C, 'local', ['-o', cur_dir + '/temp_extremal_surface.fg'])
  s = Rational(str(s.numerator)+'/'+str(s.denominator)) #this turns it into a sage thing
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
    P = Permutations(range(deg)).list()
    T = Tuples(P, rank)
    cover_rank = 1+deg*rank-deg
    num_non_connected = 0
    num_good_transfers = 0
    if verbose>1:
      print "Doing degree: ", deg
      print "There are ",T.cardinality(), " covers to do"
    if verbose > 2:
      for t in T:
        print t
        G = covering.FISubgroup(base_gens, t)
        print "Found covering subgroup ", G
        if not G.connected:
          num_non_connected += 1
          print "Not a connected cover"
          continue
        GF = F.lift(G)
        print "Found lifted fatgraph: "
        print(GF)
        GF.write_file_new('temp_lifted_fatgraph.fg')
        GFE = GF.ends()
        print "Found end list: ", GFE
        GFE_lifted = [[e.lift(G) for e in EL] for EL in GFE]
        print "Found the lifted end list: ", GFE_lifted
        compat_orders = ends.compatible_cyclic_orders(GFE_lifted, cover_rank)
        print "Found compatible cyclic orders: ", compat_orders
        if len(compat_orders) != 0:
          print "*** good transfer ", (t, G, compat_orders)
          found_transfers.append((t, G, compat_orders))
          num_good_transfers += 1
    else:
      for t in T:
        G = covering.FISubgroup(base_gens, t)
        if not G.connected:
          num_non_connected += 1
          continue
        GF = F.lift(G)
        GFE = GF.ends()
        GFE_lifted = [[e.lift(G) for e in EL] for EL in GFE]
        compat_orders = ends.compatible_cyclic_orders(GFE_lifted, cover_rank)
        if len(compat_orders) != 0:
          if verbose > 1:
            print "*** good transfer ", (t, G, compat_orders)
          found_transfers.append((t, G, compat_orders))
          num_good_transfers += 1
    if verbose>1:
      print "For this degree (", deg, "):"
      print "Num covers: ", T.cardinality()
      print "Num connected: ", T.cardinality() - num_non_connected
      print "Num transfers found: ", num_good_transfers 
  return found_transfers


















      
