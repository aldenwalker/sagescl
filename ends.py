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
    return "End(%r,%r)" % self.prefix, self.repeat
  
  def __str__(self):
    return "%r(%r)" % self.prefix, self.repeat
  
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


def find_all_suborders(EL):
  """returns a list of cyclic orders, each coming from a vertex in the 
  tree of ends"""
  #get a list of all the unique letters
  LEL = len(EL)
  gen_list = [(e[0],i) for i,e in enumerate(EL)]
  i=0
  while i<LEL-1 and gen_list[i][0] == gen_list[i+1][0]:
    i += 1
  if i == LEL-1:
    return find_all_suborders([e[1:] for e in EL] + [EL[0][0].swapcase()])
  gen_list = gen_list[i+1:] + gen_list[:i+1]
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
    new_EL = [EL[i][1::] for i in reduced_gen_dict[g]] + [g.swapcase()]
    sub_COs.append(find_all_suborders(new_EL))
  
  if [] in sub_COs:
    return []
  
  if len(our_CO) > 2:
    return [our_CO] + sub_COs
  
  return sub_COs
  
  
  
  
def compatible_cyclic_orders(EL, rank=None):
  """return a list of the cyclic orders compatible with the 
  given cyclically ordered set of ends"""
  if rank==None:
    R = word.chain_rank([e.get_word() for e in EL])
  else:
    R = rank
  
  #get all the suborders from the ends
  SO = find_all_suborders(EL)
  if SO == []:
    return []
  SO_lens = map(len, SO)
  
  #now we must assemble the suborders
  #start with the largest (why not), and then iteratively
  #try to add more by searching for other orders which contain things 
  #on the boundary
  big_ind = SO_lens.index(max(SO_lens))
  current_order = SO[ big_ind ]
  del SO[big_ind]
  while True:
    #for each pair on the boundary
    added_something = False
    for i in xrange(len(current_order)-1):
      p1 = current_order[i]
      p2 = current_order[i+1]
      #find it in the other orders
      for so in SO:
        if p1 in so and p2 in so:
          subslice = so[p1:p2][1:]
          if len(subslice) == 0:
            continue
          #test if it's ok
          for g in subslice:
            if current_order(p1, g, p2) != 0: #we can't possibly have a compatible order
              return []
          current_order.insert(p1, subslice)
          added_something=True
          break
      if added_something:
        break
    if not added_something:
      break
  
  #now we've got the largest thing we can build
  #any cyclic order which is an extension of this is acceptable
  current_order.extend_to_full_gen_order(R)
  return current_order
  
  
  
  
  