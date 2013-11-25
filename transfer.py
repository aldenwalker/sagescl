from sage.all import *
import word
import covering
import fatgraph
import ends
import scl

def rotate_min_first(L):
  LL = len(L)
  mL = min(L)
  mLi = L.index(mL)
  return L[mLi:] + L[:mLi]

def frac_to_sage_Rational(x):
  return Rational(str(x.numerator) + '/' + str(x.denominator))

def find_extremal_transfer(C_in, max_degree=None, degree_list=None, verbose=1, fatgraph_size_bound=None, just_one=False):
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
        compat_orders = ends.compatible_cyclic_orders(GFE, cover_rank)
        print "Found compatible cyclic orders: ", compat_orders
        if len(compat_orders) != 0:
          print "*** good transfer ", (t, G, compat_orders)
          GF.write_file_new('temp_good_transfer_fg.fg')
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
          if just_one:
            return found_transfers
    else:
      for t in T:
        G = covering.FISubgroup(base_gens, t)
        if not G.connected:
          num_non_connected += 1
          continue
        GF = F.lift(G)
        GFE = GF.ends()
        compat_orders = ends.compatible_cyclic_orders(GFE, cover_rank)
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
          if just_one:
            return found_transfers
    if verbose>1:
      print "For this degree (", deg, "):"
      print "Num covers: ", T.cardinality()
      print "Num connected: ", T.cardinality() - num_non_connected
      print "Num transfers found: ", num_good_transfers 
  return found_transfers


def single_transfer(C, G, verbose=1, all_orders=False):
  #get the rank
  rank, base_gens = word.chain_rank_and_gens(C)
  
  #get the cwd
  cur_dir = os.getcwd()
  
  #get an extremal surface
  s = scl.scl(C, 'local', ['-o', cur_dir + '/temp_extremal_surface.fg'])
  s = Rational(str(s.numerator)+'/'+str(s.denominator)) #this turns it into a sage thing
  F = fatgraph.read_file(cur_dir + '/temp_extremal_surface.fg')
  
  #lift the fatgraph
  GF = F.lift(G)
  if verbose>1:
    print "Found lifted fatgraph: "
    if verbose>2:
      print(GF)
    GF.write_file_new('temp_lifted_fatgraph.fg')
  GFE = GF.ends()
  if verbose>1:
    print "Found lifted end list: ", GFE
  compat_orders = ends.compatible_cyclic_orders(GFE, G.rank, all_orders=all_orders)
  if verbose>1:
    print "Found compatible cyclic orders: ", compat_orders
    GF.write_file_new('temp_good_transfer_fg.fg')
    if verbose > 1:
      print "Double checking rot = ", compat_orders[0].rot(G.lift(C))
  if len(compat_orders)>0:
    if compat_orders[0].rot(G.lift(C)) != 2*G.degree*s:
        print "Rot isn't extremal?"
        print "C = ", C
        print "G = ", G
        print "O = ", compat_orders[0]
        return []
  return compat_orders



def find_transfer_families(n, ntrials, rank=2, \
                           cover_deg_bound=4,  \
                           fatgraph_size_bound=50, \
                           covers_to_try=None,
                           verbose=1):
  found_transfer_families = []
  for i in xrange(ntrials):
    F = word.random_family(n, rank, gens=['x','y'])
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
      if cover_deg_bound > 4:
        if verbose > 1:
          print "Min cover degree ", min_cover_deg, " is too big"
        break
      
      if verbose>1:
        print "Finding extremal transfers at degree: ", [min_cover_deg]
      T = find_extremal_transfer(F(N), degree_list=[min_cover_deg], covers_to_try=covers_to_try, fatgraph_size_bound=fatgraph_size_bound, just_one=True)
      
      if len(T) > 0:
        if verbose>1:
          print "Found ", len(T), " good transfers"
        transfers.append(T)
        N += 1
      else:
        break
    
    if all([t[0][1].degree==1 for t in transfers]):
      #don't record this
      if verbose > 1:
        print "We only found extremal basic rots"
      continue
    elif len(transfers) == 0:
      if verbose>1:
        print "We found no transfers"
      continue
    
    if verbose>1:
      print "We found some good transfers"
    found_transfer_families.append( ( F, transfers ) )
      
  return found_transfer_families



def random_transfers(n, rank, ntrials, verbose=1):
  """take lots of random chains, and if the denominator is in the 
  appropriate range, try to find a transfer"""
  results_by_denominator = {}
  for i in xrange(ntrials):
    C = word.random_hom_triv_chain(n, rank)
    if verbose>1:
      print "Trying: ", C
    try:
        s = frac_to_sage_Rational(scl.scl(C))
    except:
      if verbose > 1:
        print "Scl computation failed"
      continue
    min_cover_deg = (s.denominator()/2 if s.denominator()%2 == 0 else s.denominator())
    if verbose>1:
      print "scl = ", s, "; min cover: ", min_cover_deg
    if min_cover_deg > 4:
      if verbose > 1:
        print "Min cover degree is too big"
      continue
    T = find_extremal_transfer(C, degree_list=[min_cover_deg], fatgraph_size_bound=80, just_one=True)
    if len(T)>0:
      if verbose>1:
        print "Found transfer"
      if s.denominator() in results_by_denominator:
        results_by_denominator[s.denominator()][0] += 1
        results_by_denominator[s.denominator()][1].append(C)
      else:
        results_by_denominator[s.denominator()] = [1, [C]]
    else:
      if s.denominator() in results_by_denominator:
        results_by_denominator[s.denominator()][0] += 1
      else:
        results_by_denominator[s.denominator()] = [1, []]
  return results_by_denominator


def is_cyclic_cover_order(G, O):
  """returns true iff the order corresponds to a cyclic order 
  of a bunch of copies of the rose translated along by the 
  cyclic generator"""
  if G.degree == 1:
    return False
  for i in xrange(len(G.base_gens)):
    if G.base_gen_actions[i] != range(G.degree):
      cyclic_base_gen_i = i
      break
  cyclic_base_gen = G.base_gens[cyclic_base_gen_i]
  gen_prefix_amounts = [word.count_prefixes(w, cyclic_base_gen) for w in G.gens_in_base_group()]
  gen_prefix_amounts = [gpa%G.degree for gpa in gen_prefix_amounts] #make them all positive
  order_pa_list = [gen_prefix_amounts[G.gens.index(gco.lower())] for gco in O]
  #we must ensure that the order is monotone -- it goes up, then down
  order_pa_list = rotate_min_first(order_pa_list)
  direction = 'up'
  L = len(order_pa_list)
  i=0
  while i<L-1:
    if direction=='up' and order_pa_list[i] > order_pa_list[i+1]:
      direction = 'down'
    elif direction=='down' and order_pa_list[i] < order_pa_list[i+1]:
      return False
    i+=1
  return True
      
    

def cyclic_transfer_families(n, rank, ntrials, family_bound=8, cover_degree_bound=9, verbose=1):
  gens = word.alphabet[:rank]
  found_transfers = []
  for i in xrange(ntrials):
    F = word.random_family(n, rank)
    family_transfers = []
    if verbose>1:
      print "Trying family", F
    N=1
    while True:
      try:
        s = frac_to_sage_Rational(scl.scl(F(N)))
      except:
        if verbose > 1:
          print "Scl computation failed"
          break
      min_cover_deg = (s.denominator()/2 if s.denominator()%2 == 0 else s.denominator())
      if verbose>1:
        print "N = ", N, "scl =", s, "cover degree =", min_cover_deg      
      if N>2 and min_cover_deg==1:
        if verbose>1:
          print "Family isn't getting complicated"
        break
      if N > family_bound or min_cover_deg > cover_degree_bound:
        if verbose>1:
          print "Cover is too complicated"
        break
      perms = [range(1,min_cover_deg) + [0]] + [range(min_cover_deg) for g in xrange(rank-1)]
      G = covering.FISubgroup(gens, perms)
      ET = single_transfer(F(N), G, all_orders=True)
      if verbose>1:
        print "Found all orders: ",ET
      ET = [et for et in ET if is_cyclic_cover_order(G, et)]
      if verbose>1:
        print "Found cyclic cover orders:", ET
      if len(ET) > 0:
        if verbose>1:
          print "Found transfer orders: ", ET
        family_transfers.append( (N, ET) )
        N += 1
      else:
        break
    if family_transfers != []:
      found_transfers.append( (F, family_transfers) )
  return found_transfers


def cyclic_cover_vertices(F, G):
  """returns a list of (cyclic order, path order to basepoint) for the vertices, 
  recording the cyclic order on the edges, plus the number of "y"s (or whatever 
  generator goes to the generator of the cyclic group) needed to get to the vertex"""
  F.comb()
  vert_words = [F.word_from_basepoint_to_vert(vi) for vi in xrange(len(F.V))]
  vert_orders = F.vertex_cyclic_orders()
  for i in xrange(len(G.base_gens)):
    if G.base_gen_actions[i] != range(G.degree):
      cyclic_gen = G.base_gens[i]
  vert_word_images = [(w.count(cyclic_gen) - w.count(cyclic_gen.swapcase()))%G.degree for w in vert_words]
  return [(vert_orders[i], vert_word_images[i]) for i in xrange(len(F.V))]
      
      



