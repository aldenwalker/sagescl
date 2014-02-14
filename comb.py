def comb_no_end(f_in, ntines, cmp_val='1'):
  if type(f_in) == str:
    f = open(f_in, 'w')
  else:
    f = f_in
   
  f.write( "(set-logic QF_LRA)\n")

  for i in xrange(ntines+1):
    f.write( "(declare-fun b_" + str(i) + " () Real)\n" )
  for i in xrange(2*ntines+1):
    f.write( "(declare-fun m_" + str(i) + " () Real)\n" )
  for i in xrange(2*ntines):
    f.write( "(declare-fun t_" + str(i) + " () Real)\n" )

  for i in xrange(2*ntines):
    f.write( "(declare-fun L_" + str(i) + " () Real)\n" )

  for i in xrange(ntines+1):
    f.write( "(assert (and (>= b_" + str(i) + " 0) (< b_" + str(i) + " 1)))\n" )
  for i in xrange(2*ntines+1):
    f.write( "(assert (and (>= m_" + str(i) + " 0) (< m_" + str(i) + " 1)))\n" )
  for i in xrange(2*ntines):
    f.write( "(assert (and (>= t_" + str(i) + " 0) (< t_" + str(i) + " 1)))\n" )
    
  for i in xrange(ntines):
    #bottom part
    a = "b_" + str(i)
    b = "b_" + str(i+1)
    c = "m_" + str(2*i+1)
    d = "m_" + str(2*i)
    #a<b and c<d and ((b<c and d<a+1) or (d<a and b<c+1))
    inside_or = "(or (and (< {b} {c}) (< {d} (+ 1 {a}))) (and (< {d} {a}) (< {b} (+ 1 {c}))))".format(a=a,b=b,c=c,d=d)
    clause_1 = "(and (< {a} {b}) (< {c} {d}) ".format(a=a,b=b,c=c,d=d) + inside_or + ")"
    #a<b and d<c and d<a and b<c
    clause_2 = "(and (< {a} {b}) (< {d} {c}) (< {d} {a}) (< {b} {c}))".format(a=a,b=b,c=c,d=d)
    #c<d and b<a and b<c and d<a
    clause_3 = "(and (< {c} {d}) (< {b} {a}) (< {b} {c}) (< {d} {a}))".format(a=a,b=b,c=c,d=d)
    f.write("(assert (or " + clause_1 + " " + clause_2 + " " + clause_3 + "))\n")
    
    #(a<b and L = b-a) or (b<a and L=1-(a-b)) 
    bottom = "(or (and (< {a} {b}) (= {L} (- {b} {a}))) (and (< {b} {a}) (= {L} (- 1 (- {a} {b})))))\n".format(a=a,b=b,L='L_' + str(2*i))
    #(c<d and L = c-d) or (d<c and L=1-(c-d))
    top = "(or (and (< {c} {d}) (= {L} (- {d} {c}))) (and (< {d} {c}) (= {L} (- 1 (- {c} {d})))))\n".format(c=c,d=d,L='L_' + str(2*i))
    
    f.write("(assert " + bottom + ")\n")
    f.write("(assert " + top + ")\n")
    
    #upright part
    a = "m_" + str(2*i+2)
    b = "t_" + str(2*i+1)
    c = "t_" + str(2*i)
    d = "m_" + str(2*i+1)
    #a<b and c<d and ((b<c and d<a+1) or (d<a and b<c+1))
    inside_or = "(or (and (< {b} {c}) (< {d} (+ 1 {a}))) (and (< {d} {a}) (< {b} (+ 1 {c}))))".format(a=a,b=b,c=c,d=d)
    clause_1 = "(and (< {a} {b}) (< {c} {d}) ".format(a=a,b=b,c=c,d=d) + inside_or + ")"
    #a<b and d<c and d<a and b<c
    clause_2 = "(and (< {a} {b}) (< {d} {c}) (< {d} {a}) (< {b} {c}))".format(a=a,b=b,c=c,d=d)
    #c<d and b<a and b<c and d<a
    clause_3 = "(and (< {c} {d}) (< {b} {a}) (< {b} {c}) (< {d} {a}))".format(a=a,b=b,c=c,d=d)
    f.write("(assert (or " + clause_1 + " " + clause_2 + " " + clause_3 + "))\n")
    
    #(a<b and L = b-a) or (b<a and L=1-(a-b)) 
    right = "(or (and (< {a} {b}) (= {L} (- {b} {a}))) (and (< {b} {a}) (= {L} (- 1 (- {a} {b})))))\n".format(a=a,b=b,L='L_' + str(2*i+1))
    #(c<d and L = c-d) or (d<c and L=1-(c-d))
    left = "(or (and (< {c} {d}) (= {L} (- {d} {c}))) (and (< {d} {c}) (= {L} (- 1 (- {c} {d})))))\n".format(c=c,d=d,L='L_' + str(2*i+1))
    
    f.write("(assert " + right + ")\n")
    f.write("(assert " + left + ")\n")
    
  bottom_s = '(+'
  upright_s = '(+'
  for i in xrange(ntines):
    bottom_s += ' L_' + str(2*i)
    upright_s += ' L_' + str(2*i+1)
  bottom_s += ')'
  upright_s += ')'
  
  total_s = "(+ " + upright_s + " (* 2 " + bottom_s + "))"
  
  f.write('(assert (<= ' + cmp_val + ' ' + total_s + '))\n')
  f.write('(assert (<= b_0 m_0))\n')
  f.write('(assert (<= m_0 (/ 1 10)))\n')
  f.write('(check-sat)\n')
  f.close()
  
  
    
    


def comb_smt(f_in, ntines, bound_str="1", do_proof=False):
  
  if type(f_in) == str:
    f = open(f_in, 'w')
  else:
    f = f_in
  
  
  f.write( "(set-option :produce-proofs true)\n" )
  f.write( "(set-logic QF_LRA)\n")
  
  for i in xrange(ntines+2):
    f.write( "(declare-fun b_" + str(i) + " () Real)\n" )
  for i in xrange(2*ntines+2):
    f.write( "(declare-fun m_" + str(i) + " () Real)\n" )
  for i in xrange(2*ntines):
    f.write( "(declare-fun t_" + str(i) + " () Real)\n" )

  for i in xrange(2*ntines+1):
    f.write( "(declare-fun L_" + str(i) + " () Real)\n" )

  for i in xrange(ntines+1):
    f.write( "(assert (and (>= b_" + str(i) + " 0) (< b_" + str(i) + " 1)))\n" )
  for i in xrange(2*ntines+2):
    f.write( "(assert (and (>= m_" + str(i) + " 0) (< m_" + str(i) + " 1)))\n" )
  for i in xrange(2*ntines):
    f.write( "(assert (and (>= t_" + str(i) + " 0) (< t_" + str(i) + " 1)))\n" )

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
    f.write( "(assert (or " + clause1 + " " + clause2 + "))\n" )
  
    side0 = "(or (and (<" + a + b + ") (= L_" + str(2*i) + " (- " + b + a + "))) (and (<" + b + a + ") (= L_" + str(2*i)+ " (- 1 (-" + a + b + ")))) )"
    side1 = "(or (and (<" + c + d + ") (= L_" + str(2*i) + " (- " + d + c + "))) (and (<" + d + c + ") (= L_" + str(2*i) + " (- 1 (-" + c + d + ")))) )"
    f.write( "(assert " + side0 + ")\n" )
    f.write( "(assert " + side1 + ")\n" )


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
    f.write( "(assert (or " + clause1 + " " + clause2 + "))\n" )

    side0 = "(or (and (<" + a + b + ") (= L_" + str(2*i+1) + " (- " + b + a + "))) (and (<" + b + a + ") (= L_" + str(2*i+1)+ " (- 1 (-" + a + b + ")))) )"
    side1 = "(or (and (<" + c + d + ") (= L_" + str(2*i+1) + " (- " + d + c + "))) (and (<" + d + c + ") (= L_" + str(2*i+1) + " (- 1 (-" + c + d + ")))) )"
    f.write( "(assert " + side0 + ")\n" )
    f.write( "(assert " + side1 + ")\n" )
    
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
  f.write( "(assert (or " + clause1 + " " + clause2 + "))\n" )

  side0 = "(or (and (<" + a + b + ") (= L_" + str(2*ntines) + " (- " + b + a + "))) (and (<" + b + a + ") (= L_" + str(2*ntines)+ " (- 1 (-" + a + b + ")))) )"
  side1 = "(or (and (<" + c + d + ") (= L_" + str(2*ntines) + " (- " + d + c + "))) (and (<" + d + c + ") (= L_" + str(2*ntines) + " (- 1 (-" + c + d + ")))) )"
  f.write( "(assert " + side0 + ")\n" )
  f.write( "(assert " + side1 + ")\n" )

  out_s = "(+ L_0"
  for i in xrange(ntines):
    out_s += " L_" + str(2*i+1)
  out_s += " L_" + str(2*ntines) + ")"
  in_s = "(+"
  for i in xrange(ntines-1):
    in_s += " L_" + str(2*i+2)
  in_s += ")"
  if ntines > 1:
    s = "(+ " + out_s + " (* 2 " + in_s + "))"
  else:
    s = out_s
  f.write( "(assert (<= " + bound_str + " " + s + "))\n" ) 

  f.write("(check-sat)\n")
  
  if do_proof:
    f.write( "(get-proof)\n" )
  
  f.write("(exit)\n")

  if type(f_in) == str:
    f.close()


def comb_smt_bottom(f_in, ntines, bound_str="1", do_proof=False):
  
  if type(f_in) == str:
    f = open(f_in, 'w')
  else:
    f = f_in
  
  
  f.write( "(set-option :produce-proofs true)\n" )
  f.write( "(set-logic QF_LRA)\n")
  
  for i in xrange(ntines+2):
    f.write( "(declare-fun b_" + str(i) + " () Real)\n" )
  for i in xrange(2*ntines+2):
    f.write( "(declare-fun m_" + str(i) + " () Real)\n" )

  for i in xrange(ntines+1):
    f.write( "(declare-fun L_" + str(i) + " () Real)\n" )

  for i in xrange(ntines+1):
    f.write( "(assert (and (>= b_" + str(i) + " 0) (< b_" + str(i) + " 1)))\n" )
  for i in xrange(2*ntines+2):
    f.write( "(assert (and (>= m_" + str(i) + " 0) (< m_" + str(i) + " 1)))\n" )

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
    f.write( "(assert (or " + clause1 + " " + clause2 + "))\n" )
  
    side0 = "(or (and (<" + a + b + ") (= L_" + str(i) + " (- " + b + a + "))) (and (<" + b + a + ") (= L_" + str(i)+ " (- 1 (-" + a + b + ")))) )"
    side1 = "(or (and (<" + c + d + ") (= L_" + str(i) + " (- " + d + c + "))) (and (<" + d + c + ") (= L_" + str(i) + " (- 1 (-" + c + d + ")))) )"
    f.write( "(assert " + side0 + ")\n" )
    f.write( "(assert " + side1 + ")\n" )
    
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
  f.write( "(assert (or " + clause1 + " " + clause2 + "))\n" )

  side0 = "(or (and (<" + a + b + ") (= L_" + str(ntines) + " (- " + b + a + "))) (and (<" + b + a + ") (= L_" + str(ntines)+ " (- 1 (-" + a + b + ")))) )"
  side1 = "(or (and (<" + c + d + ") (= L_" + str(ntines) + " (- " + d + c + "))) (and (<" + d + c + ") (= L_" + str(ntines) + " (- 1 (-" + c + d + ")))) )"
  f.write( "(assert " + side0 + ")\n" )
  f.write( "(assert " + side1 + ")\n" )

  s = "(* 2 (+"
  for i in xrange(ntines+1):
    s += " L_" + str(i)
  s += "))"
  f.write( "(assert (<= " + bound_str + " " + s + "))\n" ) 

  f.write("(check-sat)\n")
  
  if do_proof:
    f.write( "(get-proof)\n" )
  
  f.write("(exit)\n")

  if type(f_in) == str:
    f.close()
  










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
