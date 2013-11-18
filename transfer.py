



def find_extremal_transfer(C_in, max_degree_in=None, degree_list_in=None):
  """given a chain, try to find an extremal rot transfer.  This will 
  take a whole bunch of covers and then look through all *basic* rots 
  for all the covers"""
  
  #get a list of all the covering degrees we should go through
  if max_degree_in == None:
    if degree_list_in == None:
      print "I need a max degree or degree list"
      return None
    degree_list = degree_list_in
  else:
    degree_list = range(2,max_degree+1)
  
  #get the chain as a list of words
  if type(C_in) == str:
    C = [C_in]
  else:
    C = C_in
  
  #get the rank
  rank = word.chain_rank(C)
  
  #go through the list of degrees
  for deg in degree_list:
    #go through all possible permuations for all the generators
    #i.e. all covering spaces (ouch)
    