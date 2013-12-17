import subprocess
import random
import sys
import os
import fractions

from word import *
from HNN import min_in_orbit
import fatgraph

WDIR = './'
FNULL = open('/dev/null', 'w')
SCALLOPDIR = '/home/akwalker/Documents/software/scallop/'
TROLLOPDIR = '/home/akwalker/Documents/software/trollop/'
SCABBLEDIR = '/home/akwalker/Documents/software/scabble/'
GALLOPDIR = '/home/akwalker/Documents/software/gallop/'


def scl_surface(chain, args=[]):
  """calls scallop -local and gets the fatgraph surface"""
  
  cur_dir = os.getcwd()
  
  #get an extremal surface
  s = scl(chain, args=args + ['-o', cur_dir + '/temp_fatgraph.fg'])
  F = fatgraph.read_file(cur_dir + '/temp_fatgraph.fg')
  return F


def scl(chain, mode='local', args=[], return_fraction=True, verbose=1):

  if chain == '' or chain == ['']:
    return 0
  
  if type(chain) == str:
    C = [chain]
  else :
    C = [x for x in chain if x != '']
  
  cmd_list = ['scallop', '-' + mode] + args + C
    
  if verbose > 1:
    print "Running " + SCALLOPDIR+str(cmd_list)
  if sys.version[:3] == '2.7':
    mout = subprocess.check_output(cmd_list,  \
                             executable=SCALLOPDIR+cmd_list[0],     \
                             stderr=FNULL,                            \
                             cwd=SCALLOPDIR)
  else:
    sca = subprocess.Popen(cmd_list,  \
                       executable=SCALLOPDIR+cmd_list[0],         \
                       stdout=subprocess.PIPE,                  \
                       stderr=FNULL,                            \
                       cwd=SCALLOPDIR)
    mout = sca.communicate()[0]
  if return_fraction:
    dat = mout.split(' ')[-3]
    if dat == 'feasible':
      raise NameError('No feasible solution found')
    return fractions.Fraction(dat)
  else:
    dat = mout.split(' ')[-1]
    return eval(dat)




















########################################################################
# these functions were written with the old scallop in mind



def scylla(gen_string, C_input, method=None):
  if type(C_input) == str:
    C = [C_input]
  else :
    C = [x for x in C_input if x != '']
  
  run_string = ['scallop']
  if method != None:
    run_string += ['-m' + method]
  run_string += [gen_string]
  run_string += C
  
  if sys.version[:3] == '2.7':
    mout = subprocess.check_output(run_string,  \
                             executable=SCALLOPDIR+run_string[0],     \
                             stderr=FNULL,                            \
                             cwd=SCALLOPDIR)
  else:
    sca = subprocess.Popen(run_string,  \
                       executable=SCALLOPDIR+run_string[0],         \
                       stdout=subprocess.PIPE,                  \
                       stderr=FNULL,                            \
                       cwd=SCALLOPDIR)
    mout = sca.communicate()[0]
  dat = mout.split(' ')[-3]
  return sage_eval(dat)

def old_scl(C_input, scylla=None, scylla_i=False):
  gens =  list(set( (''.join(C_input)).lower() ))
  
  if C_input == '' or C_input == ['']:
    return 0
  
  if type(C_input) == str:
    C = [C_input]
  else :
    C = [x for x in C_input if x != '']
  
  if scylla != None:
    gen_string = ''.join( [ x + '0' for x in gens] )
    if scylla_i:
      run_string = ['scylla', '-mCIPT', gen_string] +  C
    else:
      run_string = ['scylla', gen_string] +  C
  else:
    run_string = ['scallop'] + C
    
  #print "Running " + SCALLOPDIR+run_string[0]
  if sys.version[:3] == '2.7':
    mout = subprocess.check_output(run_string,  \
                             executable=SCALLOPDIR+run_string[0],     \
                             stderr=FNULL,                            \
                             cwd=SCALLOPDIR)
  else:
    sca = subprocess.Popen(run_string,  \
                       executable=SCALLOPDIR+run_string[0],         \
                       stdout=subprocess.PIPE,                  \
                       stderr=FNULL,                            \
                       cwd=SCALLOPDIR)
    mout = sca.communicate()[0]
  dat = mout.split(' ')[-3]
  return fractions.Fraction(dat)






def scl_LP(C_input, filename, do5=False): 
  if type(C_input) == str:
    C = [C_input]
  else :
    C = [x for x in C_input if x != '']
  run_string = ['scallop']
  if do5:
    run_string += ['-m5']
  run_string += ['-L!', os.getcwd() + '/' + filename]
  run_string += C
  #print "Command: " + str(runString+['-L!', filename]+[x for x in C.split('+')])
  if sys.version[:3] == '2.7':
    mout = subprocess.check_output(run_string,  \
                             executable=SCALLOPDIR+run_string[0],     \
                             stderr=FNULL,                            \
                             cwd=SCALLOPDIR)
  else:
    sca = subprocess.Popen(run_string,  \
                       executable=SCALLOPDIR+run_string[0],         \
                       stdout=subprocess.PIPE,                  \
                       stderr=FNULL,                            \
                       cwd=SCALLOPDIR)
    mout = sca.communicate()[0]    
  return True
  
def matlab_LP(filename) :
  os.system('/usr/local/MATLAB/R2011a/bin/matlab < LP_run.m > mOut.txt')
  return True
  
def read_scl_file(filename):
  fp = open(filename, 'r')
  ft = fp.read()
  #print ft.split('\n')
  fts = ft.split('\n')
  return (float(fts[-8]), fts[-3]=='Solved')
  
  
def matlab_scl(C, do5=False) :
  scl_LP(C, 'LP', do5)
  matlab_LP('')
  data = read_scl_file('mOut.txt')
  return (data[0]/4, data[1])


def smart_scl(C_in, do5=False):
  if type(C_in) == str:
    C = [C_in]
  else:
    C = [x for x in C_in]
  C_min = min_in_orbit(C_in)
  chain_len = sum(map(len, C_min))
  if chain_rank > 2:
    if chain_len < 60:
      ans = scl( C_min, scylla=True, scylla_i=True)
    else:
      ans = None
  else:
    if chain_len < 56:
      ans = scl( C_min )
    elif chain_len < 90:
      ans = matlab_scl( C_min )[0]
    else:
      ans = None
  return ans


def trollop(C_input, ell, rank=None, scl=False, mat_comp=False, sep_domain=False, method='GLPK'):
  if not mat_comp:
    if type(C_input) == str:
      C = [C_input]
    else:
      C = [x for x in C_input if x != '']
    run_string = ['trollop']
    if scl:
      run_string += ['-scl']
    run_string += [str(ell)]
    run_string += C
  else:
    run_string = ['trollop']
    if method != 'GLPK':
      run_string.append('-M'+method)
    if sep_domain:
      run_string.append('-dom')
    run_string.extend(['-mat', 'M.mat', 'N.mat', 'b.vec', str(rank), str(ell)]) 
  
  
  if sys.version[:3] == '2.7':
    mout = subprocess.check_output(run_string,  \
                                   executable=TROLLOPDIR+run_string[0],     \
                                   stderr=FNULL,                            \
                                   cwd=TROLLOPDIR)
  else:
    sca = subprocess.Popen(run_string,  \
                       executable=TROLLOPDIR+run_string[0],         \
                       stdout=subprocess.PIPE,                  \
                       stderr=FNULL,                            \
                       cwd=TROLLOPDIR)
    mout = sca.communicate()[0]      
  dat = mout.split()[-3]
  return fractions.Fraction(dat)
  

def smallest_positive_scl_vertex(C1_in, C2_in):
  if type(C1_in) == str:
    C1 = [C1_in]
  else:
    C1 = C1_in
  if type(C2_in) == str:
    C2 = [C2_in]
  else:
    C2 = C2_in
  
  run_string = ['scabble', '-L+', '/dev/null']
  run_string += C1
  run_string += [',']
  run_string += C2
  if sys.version[:3] == '2.7':
    mout = subprocess.check_output(run_string,  \
                                   executable=SCABBLEDIR+run_string[0],     \
                                   stderr=FNULL,                            \
                                   cwd=SCABBLEDIR)
  else:
    sca = subprocess.Popen(run_string,  \
                       executable=SCABBLEDIR+run_string[0],         \
                       stdout=subprocess.PIPE,                  \
                       stderr=FNULL,                            \
                       cwd=SCABBLEDIR)
    mout = sca.communicate()[0]      
  dat = mout.split(',')
  return map(fractions.Fraction, dat)
  
def extremal_surface(chain, weights_in=None):
  if weights_in == None:
    weights = [1 for x in chain]
  else:
    weights = weights_in
    
  run_string = ['scallop','-s', os.getcwd() + '/surf']
  for i in xrange(len(chain)):
    run_string += [str(weights[i])+chain[i]]
      
  if sys.version[:3] == '2.7':
    mout = subprocess.check_output(run_string,  \
                             executable=SCALLOPDIR+run_string[0],     \
                             stderr=FNULL,                            \
                             cwd=SCALLOPDIR)
  else:
    sca = subprocess.Popen(run_string,  \
                       executable=SCALLOPDIR+run_string[0],         \
                       stdout=subprocess.PIPE,                  \
                       stderr=FNULL,                            \
                       cwd=SCALLOPDIR)
    mout = sca.communicate()[0]
  
  #read the extremal surface
  return fatgraph.read_file(os.getcwd() + '/surf.fg')
  
  
  
def gallop(C_in, 
           folded=False, 
           ffolded=None, 
           only_check_exists=False,
           trivalent=False,
           pSides = None,
           save_file=None,
           time_limit=0,
           solver="GLPK"):
  if type(C_in) == str:
    C = [C_in]
  else:
    C = C_in
  #if type(weights_in) == int:
  #  weights = [weights_in]
  #if weights != None:
  #  C = [str(weights[i]) + C[i] for i in xrange(len(C))]
  run_string = ['scallop', '-local']
  if folded:
    run_string.append('-f')
  if ffolded != None:
    if ffolded == True:
      run_string.append('-ff')
    else:
      run_string.append('-ff' + str(ffolded))
  if only_check_exists:
    run_string.append('-e')
  if pSides != None:
    run_string.append('-p' + str(pSides))
  if trivalent:
    run_string.append('-t')
  if solver != 'GLPK':
    if time_limit != 0:
      run_string.append('-mGUROBI')
      run_string.append('-t'+ str(time_limit))
    else:
      run_string.append('-mGUROBI')
  if save_file != None:
    run_string.append('-o')
    run_string.append(save_file)
  run_string += C
  if sys.version[:3] == '2.7':
    mout = subprocess.check_output(run_string,  \
                                   executable=SCALLOPDIR+run_string[0],     \
                                   stderr=FNULL,                            \
                                   cwd=SCALLOPDIR)
  else:
    sca = subprocess.Popen(run_string,  \
                       executable=SCALLOPDIR+run_string[0],         \
                       stdout=subprocess.PIPE,                  \
                       stderr=FNULL,                            \
                       cwd=SCALLOPDIR)
    mout = sca.communicate()[0]      
  dat = mout.split('\n')
  if only_check_exists:
    if dat[0].split(' ')[0] == "Feasible":
      return True
    else:
      return False
  return fractions.Fraction(dat[0].split(' ')[-3])
  












