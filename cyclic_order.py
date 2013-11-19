import word

class CyclicOrder:
  def __init__(self, inp):
    if type(inp) == str:
      self.CO = inp
    else:
      self.CO = ''.join(inp)
    self.L = len(self.CO)
    self.indices = dict([(g,i) for (i,g) in enumerate(self.CO)])
    
  
  def __repr__(self):
    return "CyclicOrder(" + self.CO + ")"
  
  def __str__(self):
    return self.CO
  
  def __getitem__(self, n):
    if type(n) == slice:
      if type(n.start) == str:
        start = self.indices[n.start]
        stop = self.indices[n.stop]
      else:
        start = n.start
        stop = n.stop
      step = n.step
      if stop < start:
        return self.CO[start::step] + self.CO[:stop:step]
      else:
        return self.CO[start:stop:step]
    return self.CO[n%self.L]
  
  def __contains__(self, x):
    return x in self.indices
  
  def insert(self, location, subslice):
    if type(location) == str:
      loc = self.indices[location]
    else:
      loc = location
    self.__init__(self.CO[:loc+1] + subslice + self.CO[loc+1:])
   
  def __call__(self, a, b, c):
    ai = self.indices.get(a, None)
    if ai == None:
      return 0
    bi = self.indices.get(b, None)
    if bi == None:
      return 0
    ci = self.indices.get(c, None)
    if ci == None:
      return 0
    if bi < ai:
      bi += self.L
    if ci < ai:
      ci += self.L
    return (1 if bi < ci else -1)
  
  def __len__(self):
    return self.L
  
  def extend_to_full_gen_order(self, rank):
    """extend self to a full order on all generators"""
    gens = word.alphabet[:rank]
    gens += [x.swapcase() for x in gens]
    gens_to_add = [g for g in gens if g not in self]
    gens_to_add = ''.join(gens_to_add)
    self.__init__(self.CO + gens_to_add)
    