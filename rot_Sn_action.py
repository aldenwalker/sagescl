from word import *



class CQ:
  def __init__(self, rules=None):
    if rules == None:
      self.rules = {}
      self.lens = []
    else:
      self.rules = rules
      self.lens = list(set(map(len,setself.rules.keys())))
      self.cleanup()

  def ap(self, w, kind='chain'):
    if type(w) == list:
      return sum([self.ap(x, kind) for x in w])
    Lw = len(w)
    s = 0
    if kind == 'chain':
      for i in xrange(Lw):
        for j in self.lens:
          x = cyclic_subword(w, i, j)
          s += self.rules.get(x, 0)
      return s
    elif kind == 'word':
      for i in xrange(Lw):
        for j in self.lens:
          if i+j > Lw:
            continue
          x = cyclic_subword(w, i, j)
          s += self.rules.get(x, 0)
      return s
    else:
      print "kind not recognized"
  
  def add_words(self, d):
    




#a permutation class
class Perm:
  identity_rules = string.maketrans('','')
  
  def __init__(self, rules=['',''], extend_random=True):
    if type(rules) == list:
      targets = [x for x in letters if x not in rules[1]]
      sources = [x for x in letters if x not in rules[0]]
      if extend_random:
        random.shuffle(targets)
      self.rules = string.maketrans(rules[0] + ''.join(sources), \
                                    rules[1] + ''.join(targets))
    elif type(rules) == str:
      self.rules = rules
  
  def __repr__(self):
    return self.rules[97:123]+self.rules[32]
  
  def __str__(self):
    ans = ''.join([a + ' ' for a in alph])
    ans += '\n' + (27*'| ') + '\n'
    targets = self.rules[97:123] + self.rules[32]
    ans += ''.join([x+' ' for x in targets])
    return ans
  
  def copy(self):
    return Perm(self.rules)
  
  def encode(self, s):
    return string.translate(s, self.rules)
  
  def inverse(self):
    target = string.maketrans('','')
    new_rules = string.maketrans(self.rules, target)
    return Perm(new_rules)
  
  #swap the letters i and j; note space is letter 26
  def swap_letters(self, i, j):
    if i == 26:
      i = 32
    else:
      i += 97
    if j == 26:
      j = 32
    else:
      j += 97
    i, j = min(i,j), max(i,j)
    self.rules = ''.join([self.rules[:i],    \
                          self.rules[j],     \
                          self.rules[i+1:j], \
                          self.rules[i],     \
                          self.rules[j+1:]])
