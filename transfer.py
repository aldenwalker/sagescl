from sage.all import *
import word
import covering
import fatgraph
import ends
import scl

def frac_to_sage_Rational(x):
  return Rational(str(x.numerator) + '/' + str(x.denominator))

def find_extremal_transfer(C_in, max_degree=None, degree_list=None, verbose=1, fatgraph_size_bound=None):
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
  
  if fatgraph_size_bound != None and len(F.V) > fatgraph_size_bound:
    if verbose>1:
      print "The original fatgraph is above the size bound"
    return []

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
          if verbose > 1:
            print "Double checking rot = ", compat_orders[0].rot(G.lift(C))
          if compat_orders[0].rot(G.lift(C)) != 2*deg*s:
            print "Rot isn't extremal?"
            print "C = ", C
            print "G = ", G
            print "O = ", compat_orders[0]
            return []
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
          found_transfers.append((t, G, compat_orders))
          if verbose > 1:
            print "*** good transfer ", (t, G, compat_orders)
            print "Double checking rot = ", compat_orders[0].rot(G.lift(C))
          if compat_orders[0].rot(G.lift(C)) != 2*deg*s:
            print "Rot isn't extremal?"
            print "C = ", C
            print "G = ", G
            print "O = ", compat_orders[0]
            return []
          num_good_transfers += 1
    if verbose>1:
      print "For this degree (", deg, "):"
      print "Num covers: ", T.cardinality()
      print "Num connected: ", T.cardinality() - num_non_connected
      print "Num transfers found: ", num_good_transfers 
  return found_transfers


def find_transfer_families(n, ntrials, rank=2, verbose=1):
  found_transfer_families = []
  for i in xrange(ntrials):
    F = word.random_family(n, rank)
    if verbose>1:
      print "\nTrying family: ", F
    transfers = []
    N = 1
    while True:
      try:
        s = frac_to_sage_Rational(scl.scl(F(N)))
      except:
        if verbose > 1:
          print "Scl computation failed"
        break
      min_cover_deg = (s.denominator()/2 if s.denominator()%2 == 0 else s.denominator())
      if verbose > 1:
        print "Trying F(" +  str(N) +  ") = " + str(F(N))
        print "scl = ", s
        print "Min cover degree: ", min_cover_deg
      
      if N>2 and min_cover_deg == 1:
        if verbose > 1:
          print "Family scl isn't getting complicated"
        break
      if N>10:
        if verbose>1:
          print "Family is getting too complicated"
        break
      if min_cover_deg > 4:
        if verbose > 1:
          print "Min cover degree ", min_cover_deg, " is too big"
        break
      
      if verbose>1:
        print "Finding extremal transfers at degree: ", [min_cover_deg]
      T = find_extremal_transfer(F(N), degree_list=[min_cover_deg], fatgraph_size_bound=50)
      
      if len(T) > 0:
        if verbose>1:
          print "Found ", len(T), " good transfers"
        transfers.append(T)
        N += 1
      else:
        break
    
    if len(transfers) == 1 and transfers[0][0][1].degree == 1:
      #don't record this
      if verbose > 1:
        print "We only found an extremal basic rot"
      continue
    elif len(transfers) == 0:
      if verbose>1:
        print "We found no transfers"
      continue
    
    if verbose>1:
      print "We found some good transfers"
    found_transfer_families.append( ( F, transfers ) )
      
  return found_transfer_families
















      
