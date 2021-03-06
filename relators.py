from word import *

def reduce_disk( D ) :
  """reduce a disk -- i.e. reduce the boundary by folding"""
  LD = [D[0], list(D[1])]
  while True:
    nri = non_reduced_index( LD[0] )
    if nri == None:
      break
    #just remove the letters at nri and nri+1
    LLD0 = len(LD[0])
    if nri == LLD0-1:
      LD[0] = LD[0][1:LLD0-1]
      LD[1] = LD[1][1:LLD0-1]
    elif nri == 0:
      LD[0] = LD[0][2:]
      LD[1] = LD[1][2:]
    elif nri == LLD0-2:
      LD[0] = LD[0][:LLD0-2]
      LD[1] = LD[1][:LLD0-2]
    else:
      LD[0] = LD[0][:nri] + LD[0][ (nri+2): ]
      LD[1] = LD[1][:nri] + LD[1][ (nri+2): ]
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

def all_ways_to_add_relator( (db, lets), R ):
  """returns a set of disks which give all outcomes from 
  adding one of the relators to D = (d,lets).  It won't add a relator to 
  its inverse in the matching location.  The relators must be listed [relators, inverses]"""
  new_disks = set()
  #print "Adding to the disk ", (db, lets)
  for ri, r in enumerate(R): #try every relator
    LR = len(R)
    Lr = len(r)
    lets_in_r = letter_index_dict(r)
    relator_disk = ( r, tuple([ (ri, j) for j in xrange(len(r))]) ) #the disk which is the relator
    rii = (ri + (LR/2)) % LR  # the index of the inverse relator
    #print "Adding relator: ", relator_disk
    #print "With letter indices: ", lets_in_r
    #print "And inverse relator ", rii
    for dbi in xrange(len(db)):
      #print "Trying spot ", dbi
      relator_gluing_indices = lets_in_r.get( inverse(db[dbi]), [] )
      #print "It can glue on at indices ", relator_gluing_indices
      for rgi in relator_gluing_indices:
        #print "Thinking about gluing new relator letter ", relator_disk[1][rgi], " to ", lets[dbi]
        if lets[dbi][0] == rii and lets[dbi][1] == (Lr - 1) - rgi:
          #print "Nope this would be nonreduced"
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
    lcr = least_cyclic_rep(r)
    if lcr in R or inverse(lcr) in R:
      continue
    R.append(lcr)
  RI = [inverse(r) for r in R]
  R = R + RI

  old_disks = set()
  current_disks = set([ (r, tuple([(i,j) for j in xrange(len(r))])) for i,r in enumerate(R)])
  current_disk_size = 1
  while current_disk_size < max_num_tiles:
    new_disks = set()
    for D in current_disks:
      new_disks = new_disks.union( all_ways_to_add_relator( D, R ) )
    #print "Current disk size: ", current_disk_size
    #print "Had the current disks: "
    #print current_disks
    #print "Just made the new disks: "
    #print new_disks
    old_disks = old_disks.union(current_disks)
    current_disks = new_disks
    current_disk_size += 1
  return [ r for (r,lets) in (old_disks | current_disks)]
