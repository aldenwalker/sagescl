

def comb_smt(ntines):
  for i in xrange(ntines+2):
    print "(declare-fun b_" + str(i) + " () Real)"
  for i in xrange(2*ntines+2):
    print "(declare-fun m_" + str(i) + " () Real)"
  for i in xrange(2*ntines):
    print "(declare-fun t_" + str(i) + " () Real)"

  for i in xrange(2*ntines+1):
    print "(declare-fun L_" + str(i) + " () Real)"

  for i in xrange(ntines+1):
    print "(assert (and (>= b_" + str(i) + " 0) (< b_" + str(i) + " 1)))"
  for i in xrange(2*ntines+2):
    print "(assert (and (>= m_" + str(i) + " 0) (< m_" + str(i) + " 1)))"
  for i in xrange(2*ntines):
    print "(assert (and (>= t_" + str(i) + " 0) (< t_" + str(i) + " 1)))"

  #do overlaps in each L shape for each tine
  for i in xrange(ntines):
    #bottom of the L
    a = " b_" + str(i) + " "
    b = " b_" + str(i+1) + " "
    c = " m_" + str(2*i+1) + " "
    d = " m_" + str(2*i) + " "
    #a<b and c<d and ((b<c and d<a+1) or (d<a and b<c+1))
    small_or = "(or (and (<" + b + c + ") (<" + d + "(+ 1" + a + "))) (and (<" + d + a + ") (<" + b + "(+ 1" + c + "))))"
    clause1 = "(and (<" + a + b + ") (<" + c + d + ") " + small_or + ")"
    #a<b and d<c and d<a and b<c
    clause2 = "(and (<" + a + b + ") (<" + d + c + ") (<" + d + a + ") (<" + b + c + "))"
    print "(assert (or " + clause1 + " " + clause2 + "))"
  
    side0 = "(or (and (<" + a + b + ") (= L_" + str(2*i) + " (- " + b + a + "))) (and (<" + b + a + ") (= L_" + str(2*i)+ " (- 1 (-" + a + b + ")))) )"
    side1 = "(or (and (<" + c + d + ") (= L_" + str(2*i) + " (- " + d + c + "))) (and (<" + d + c + ") (= L_" + str(2*i) + " (- 1 (-" + c + d + ")))) )"
    print "(assert " + side0 + ")"
    print "(assert " + side1 + ")"


    #upright part of the L
    a = " m_" + str(2*(i+1)) + " "
    b = " t_" + str(2*i+1) + " "
    c = " t_" + str(2*i) + " "
    d = " m_" + str(2*i+1) + " "
    #a<b and c<d and ((b<c and d<a+1) or (d<a and b<c+1))
    small_or = "(or (and (<" + b + c + ") (<" + d + "(+ 1" + a + "))) (and (<" + d + a + ") (<" + b + "(+ 1" + c + "))))"
    clause1 = "(and (<" + a + b + ") (<" + c + d + ") " + small_or + ")"
    #a<b and d<c and d<a and b<c
    clause2 = "(and (<" + a + b + ") (<" + d + c + ") (<" + d + a + ") (<" + b + c + "))"
    print "(assert (or " + clause1 + " " + clause2 + "))"

    side0 = "(or (and (<" + a + b + ") (= L_" + str(2*i+1) + " (- " + b + a + "))) (and (<" + b + a + ") (= L_" + str(2*i+1)+ " (- 1 (-" + a + b + ")))) )"
    side1 = "(or (and (<" + c + d + ") (= L_" + str(2*i+1) + " (- " + d + c + "))) (and (<" + d + c + ") (= L_" + str(2*i+1) + " (- 1 (-" + c + d + ")))) )"
    print "(assert " + side0 + ")"
    print "(assert " + side1 + ")"
    
  #the last edge
  a = " b_" + str(ntines) + " " 
  b = " b_" + str(ntines+1) + " "
  c = " m_" + str(2*ntines+1) + " "
  d = " m_" + str(2*ntines) + " "
  #a<b and c<d and ((b<c and d<a+1) or (d<a and b<c+1))
  small_or = "(or (and (<" + b + c + ") (<" + d + "(+ 1" + a + "))) (and (<" + d + a + ") (<" + b + "(+ 1" + c + "))))"
  clause1 = "(and (<" + a + b + ") (<" + c + d + ") " + small_or + ")"
  #a<b and d<c and d<a and b<c
  clause2 = "(and (<" + a + b + ") (<" + d + c + ") (<" + d + a + ") (<" + b + c + "))"
  print "(assert (or " + clause1 + " " + clause2 + "))"

  side0 = "(or (and (<" + a + b + ") (= L_" + str(2*ntines) + " (- " + b + a + "))) (and (<" + b + a + ") (= L_" + str(2*ntines)+ " (- 1 (-" + a + b + ")))) )"
  side1 = "(or (and (<" + c + d + ") (= L_" + str(2*ntines) + " (- " + d + c + "))) (and (<" + d + c + ") (= L_" + str(2*ntines) + " (- 1 (-" + c + d + ")))) )"
  print "(assert " + side0 + ")"
  print "(assert " + side1 + ")"

  out_s = "(+ L_0"
  for i in xrange(ntines):
    out_s += " L_" + str(2*i+1)
  out_s += " L_" + str(2*ntines) + ")"
  in_s = "(+"
  for i in xrange(ntines-1):
    in_s += " L_" + str(2*i+2)
  in_s += ")"
  s = "(+ " + out_s + " (* 2 " + in_s + "))"
  print "(assert (<= 1 " + s + "))" 


  


  










def comb_SMT_old(ntines):
  for i in xrange(ntines+2):
    print "(declare-fun p_" + str(i) + " () Real)"
  for i in xrange(2*ntines + 2):
    print "(declare-fun ell_" + str(i) + " () Real)"

  for i in xrange(ntines+2):
    print "(assert (>= p_" + str(i) + " 0))"
    print "(assert (< p_" + str(i) + " 1))"

  for i in xrange(2*ntines + 2):
    print "(assert (>= ell_" + str(i) + " 0))"

  #the vertical edges
  #[p_i,                  p_i + ell_{2i+1}]
  #[p_{i+1} + ell_{2i+3} + ell_{2i+2}, p_{i+1} + ell_{2i+3} + ell_{2i+2} + ell_{2i+1}]
  #for 0<=i<ntines (forget ell_2i+3 if i==ntines-1
  for i in xrange(ntines):
    a = "p_" + str(i)
    b = "(+ p_" + str(i) + " ell_" + str(2*i+1) + ")"
    if i<ntines-1:
      c = "(+ p_" + str(i+1) + " (+ ell_" + str(2*i+3) + " ell_" + str(2*i+2) + "))"
      d = "(+ p_" + str(i+1) + " (+ ell_" + str(2*i+3) + " (+ ell_" + str(2*i+2) + " ell_" + str(2*i+1) + ")))"
    else:
      c = "(+ p_" + str(i+1) + " ell_" + str(2*i+2) + ")"
      d = "(+ p_" + str(i+1) + " (+ ell_" + str(2*i+2) + " ell_" + str(2*i+1) + "))"
  #could have d<a
  option0 = "(and (< " + d + " " + a + ") \n   (< " + b + " " + "(+ 1 " + c + ")))"
  #could have b<c and d < a+1 
  option1 = "(and (< " + b + " " + c + ") \n   (< " + d + " " + "(+ 1 " + a + ")))"
  #could have "b+1 < c and d < a+2
  option2 = "(and (< (+ 1 " + b + ") " + c + ") \n    (< " + d + " " + "(+ 2 " + a + ")))"
  print "(or \n  (or \n   " + option0 + "\n   " + option1 + ")\n   " + option2 +")"

  #now do the bottom tines
  last_i = ntines+1
  #[p_last_i + sum_{j=0}^i ell_{2j}, p_last_i + sum_{j=0}^
