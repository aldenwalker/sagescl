def non_reduced_index( w ):
  LW = len(w)
  for i in xrange(LW-1):
    if w[i] == w[i+1].swapcase():
      return i
  if w[LW-1] == w[0].swapcase():
    return LW-1
  return None

def reduce_disk( D ) :
  """reduce a disk -- i.e. reduce the boundary by folding"""
  LD = [D[0], list(D[1])]
  while True:
    nri = non_reduced_index( LD[0] )
    if nri == None:
      break
    #just remove the letters at nri and nri+1
    LD[0] = LD[0][:nri] + LD[0][ ((nri+1)%len(LD[0])): ]
    LD[1] = LD[1][:nri] + LD[1][ ((nri+1)%len(LD[1])): ]
  return (LD[0], tuple(LD[1]))

def glue_disks_together((d1_word, d1_lets), d1_ind, (d2_word,d2_lets), d2_ind):
  """glue the disk (d1_word, d1_lets) to the disk (d2_word, d2_lets) along the 
  letter at index d1_ind in D1 and at index d2_ind in D2.  
  It folds all adjacent letters to make the result reduced"""
  #we want to remove the letters at d1_ind and d2_ind, patch the lists together, 
  #and then fold 
  glued_word = d1_word[:d1_ind] + d2_word[(d2_ind+1):] + d2_word[:d2_ind] + d1_word[(d1_ind+1):]
  glued_lets = d1_lets[:d1_ind] + d2_lets[(d2_ind+1):] + d2_lets[:d2_ind] + d1_lets[(d1_ind+1):]
  return reduce_disk( (glued_word, glued_lets) )

def letter_index_dict(w):
  ans = {}
  for let, i in enumerate(w):
    ans[let] = ans.get(let, []) + [i]
  return ans

def all_ways_to_add_relator( (db, lets), R ):
  """returns a set of disks which give all outcomes from 
  adding one of the relators to D = (d,lets).  It won't add a relator to 
  its inverse in the matching location.  The relators must be listed [relators, inverses]"""
  new_disks = set()
  for r,ri in enumerate(R): #try every relator
    LR = len(R)
    lets_in_r = letter_index_dict(r)
    relator_disk = ( r, tuple([ (ri, j) for j in xrange(LR)]) ) #the disk which is the relator
    rii = (ri + (LR/2)) % len(LR)  # the index of the inverse relator
    for dbi in xrange(len(db)):
      relator_gluing_indices = lets_in_r[ inverse(db[dbi]) ]
      for rgi in relator_gluing_indices:
        if lets[dbi][0] == rii and lets[dbi][1] == LR - rgi:
          continue  #this would be a non-reduced van-Kampen diagram
        new_disks.add( glue_disks_together( relator_disk, rgi, (db,lets), dbi ) )
  
  #now new_disks should have all disks which are one-larger than (db,lets)
  return new_disks


def gen_disk_boundaries(R_in, max_num_tiles):
  """generate all disk boundaries (trivial words) in the relators given by R"""
  #a disk is a pair (b, [(i_1, j_1), (i_2, j_2), ...]), where b is the boundary, 
  #and (i,j) records which relator (i) and letter (j).  This way, we never 
  #glue a relator to an inverse relator
  R = []
  for r in R_in:
    mcr = min_cyclic(r)
    if mcr in R or inverse(mcr) in R:
      continue
    R.append(mcr)
  RI = [inverse(r) for r in R]
  R = R + RI

  old_disks = []
  current_disks = set([ (r, tuple([(i,j) for j in xrange(len(r))])) for r,i in enumerate(R)])
  current_disk_size = 1
  while current_disk_size < max_num_tiles:
    new_disks = []
    for D in current_disks:
      new_disks = new_disks.union( all_ways_to_add_relator( D, R ) )
    old_disks = old_disks.union(current_disks)
    current_disks = new_disks
  return [ r for (r,lets) in (old_disks | current_disks)]
