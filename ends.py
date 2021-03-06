import sage.all as SAGE

import word
import cyclic_order

#an end is a finite prefix, followed by an infinitely repeating word
class FreeGroupEnd:
  def __init__(self, prefix, repeat):
    self.prefix = prefix
    self.pl = len(self.prefix)
    self.repeat = repeat
    self.rl = len(self.repeat)
  
  def __repr__(self):
    return "End(%r,%r)" % (self.prefix, self.repeat)
  
  def __str__(self):
    return "%r(%r)" % (self.prefix, self.repeat)
  
  def get_word(self, n=1):
    return self.prefix + (n*self.repeat)
  
  def get_letter(self, n):
    if n < self.pl:
      return self.prefix[n]
    else:
      return self.repeat[(n-self.pl)%self.rl]
  
  def __getitem__(self, n):
    if type(n) == slice:
      if n.stop == None:
        if n.start >= self.pl:
          if (n.start-self.pl)%self.rl != 0:
            return FreeGroupEnd(self.repeat[(n.start-self.pl)%self.rl:], self.repeat)
          else:
            return FreeGroupEnd('', self.repeat)
        else:
          return FreeGroupEnd(self.prefix[n.start:], self.repeat)
      else:
        if n.step != None:
          return ''.join([self.get_letter(i) for i in xrange(n.start, n.stop, n.step)])
        else:
          return ''.join([self.get_letter(i) for i in xrange(n.start, n.stop)])
    return self.get_letter(n)
  
  def simplify(self):
    self.prefix = word.word_reduce(self.prefix)
    self.repeat = word.word_reduce(self.repeat)
    core, core_prefix = word.cyc_red_get_conjugate(self.repeat)
    ppre = word.multiply_words(self.prefix, core_prefix)
    if ppre == '':
      #print "Simplified ", str(self.prefix), str(self.repeat)
      #print "To ", ppre, core
      return FreeGroupEnd('', core)
    while len(ppre) > 0 and ppre[-1] == core[0].swapcase():
      ppre = word.multiply_words(ppre, core)
    #print "Simplified ", str(self.prefix), str(self.repeat)
    #print "To ", ppre, core
    return FreeGroupEnd(ppre, core)

  def lift(self, G):
    """lift the end to the FISubgroup G"""
    new_pref, dest_vert = G.rewrite_path_from_base_gens_get_end(self.prefix) 
    new_rep = G.rewrite_path_from_base_gens(self.repeat, dest_vert)
    return FreeGroupEnd(new_pref, new_rep).simplify()

  def ap_morph(self, phi):
    ppre = phi.ap(self.prefix)
    prea = phi.ap(self.repeat)
    core, core_prefix = word.cyc_red_get_conjugate(prea)
    ppre = multiply_words(ppre, core_prefix)
    if ppre == '':
      return FreeGroupEnd('', core)
    while ppre[-1] == core[0].swapcase():
      ppre += core
    return FreeGroupEnd(ppre, core)


def find_all_suborders(EL):
  """returns a list of cyclic orders, each coming from a vertex in the 
  tree of ends"""
  
  #first, if all the first letters are the same, remove them
  LEL = len(EL)
  if LEL == 1:
    return []
  i=0
  while all([EL[0][i]==e[i] for e in EL]):
    i += 1
  #print "Found they agree for ", i, " characters"
  if i > 0:
    return find_all_suborders([e[i:None] for e in EL] + [EL[0][i-1].swapcase()])
  
  #now shift the ends around so that no letter shows up at the end and 
  #beginning
  gen_list = [(e[0],i) for i,e in enumerate(EL)]
  i=0
  while i<LEL-1 and gen_list[i][0] == gen_list[i+1][0]:
    i += 1
  gen_list = gen_list[i+1:] + gen_list[:i+1]
  
  #print "Gen list: ", gen_list
  
  #now the gen list ends and begins with different gens, so each 
  #should show up in at most one run
  reduced_gen_list = []
  reduced_gen_dict = dict()
  for (g,i) in gen_list:
    if g in reduced_gen_dict:
      if reduced_gen_dict[g][-1] != i-1:
        return []
      reduced_gen_dict[g].append(i)
    else:
      reduced_gen_list.append(g)
      reduced_gen_dict[g] = [i]
  #first, we have our suborder, given by reduced_gen_list
  our_CO = cyclic_order.CyclicOrder(reduced_gen_list)
  
  #the sub COs come from each of the subtrees
  sub_COs = []
  for g in reduced_gen_dict:
    if len(reduced_gen_dict[g]) == 1:
      continue
    new_EL = [EL[i] for i in reduced_gen_dict[g]]
    sub_CO = find_all_suborders(new_EL)
    if sub_CO == []:
      return []
    sub_COs += sub_CO
  
  if len(our_CO) > 2:
    return [our_CO] + sub_COs
  
  return sub_COs
  
  
  
  
def compatible_cyclic_orders(EL, rank=None, all_orders=False):
  """return a list of the cyclic orders compatible with the 
  given cyclically ordered set of ends"""
  if rank==None:
    if type(EL[0]) == list:
      R = max([word.chain_rank([e.get_word() for e in E]) for E in EL])
    else:
      R = word.chain_rank([e.get_word() for e in EL])
  else:
    R = rank
  
  #get all the suborders from the ends
  #we accept a list of ends or a list of lists of ends
  if type(EL[0]) == list:
    SO = []
    for e in EL:
      local_orders = find_all_suborders(e)
      if local_orders == []:
        return []
      SO += local_orders
  else:
    SO = find_all_suborders(EL)
  if SO == []:
    return []
  
  #print "Found all the suborders: ", SO
  
  #we check all the 4-tuples
  gens = word.alphabet[:R]
  gens += [g.swapcase() for g in gens]
  unknown_4_subsets = list(SAGE.Subsets(gens, 4))
  known_4_tuple_orders = []
  
  #get all the 4-tuples that we can from the orders that we have
  for i in xrange(len(unknown_4_subsets)-1, -1, -1):
    #print "I'm checking if we know the order on: ", unknown_4_subsets[i]
    got_order = cyclic_order.four_tuple_from_cyclic_orders(unknown_4_subsets[i], SO)
    #print "The order is: ", str(got_order)
    if got_order == None:
      return []
    if got_order != False:
      del unknown_4_subsets[i]
      known_4_tuple_orders.append(got_order)
  
  #now we need to build up the consequences
  #by continuing to scan through the unknown guys until we can conclude nothing else
  while len(known_4_tuple_orders)>0:
    did_something = False
    for i in xrange(len(unknown_4_subsets)-1, -1, -1):
      got_order = cyclic_order.four_tuple_from_cyclic_orders(unknown_4_subsets[i], known_4_tuple_orders)
      if got_order == None:
        return []
      if got_order != False:
        did_something = True
        del unknown_4_subsets[i]
        known_4_tuple_orders.append(got_order)
    if not did_something:
      break

  #print "Found the known 4 tuples: ", known_4_tuple_orders
  
  #determine if the order can be extended (any contradictions)
  tripods = list(SAGE.Subsets(gens, 3))
  SO_tripods = [O for O in SO if len(O)==3]
  everything_known = SO_tripods + known_4_tuple_orders
  for t in tripods:
    if None == cyclic_order.multiple_cyclic_order_eval(t, everything_known):
      return []

  #now we need to extend what we have to a cyclic order
  #there might be some tripods in the original thing, 
  #and we should include those
  if all_orders:
    CO = cyclic_order.all_order_extensions(R, SO_tripods + known_4_tuple_orders)
  else:
    CO = [ cyclic_order.extend_suborders_to_order(R, SO_tripods + known_4_tuple_orders) ]
  #CO = list(set(CO))
  return CO
  
  
  
  
  
  
  
def plot_free_group(rank, paths=None, edges=None, ends=None, ends_with_basepoints=None):
  G = SAGE.Graphics()
  G += SAGE.circle((0,0),1)
  G.show()
  return
  
  
  
  
  
