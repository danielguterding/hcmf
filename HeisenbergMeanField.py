import os
import subprocess
import numpy as np
import scipy.sparse as spsparse
import scipy.sparse.linalg as splinalg
import itertools

class SpinState:
  def __init__(self, statevec):
    self.statevec = np.array([bool(int(x)) for x in statevec], dtype=bool)
  def get_sz_site_resolved(self):
    return self.statevec-0.5*np.ones(self.statevec.shape)

class HeisenbergBond:
  def __init__(self, s1, s2, Jidx):
    self.s1 = int(s1)
    self.s2 = int(s2)
    self.Jidx = int(Jidx)

class MeanFieldBond:
  def __init__(self, s1, s2, Jidx):
    self.sr = int(s1) #real site 
    self.sm = int(s2) #mirror site outside of cluster
    self.Jidx = int(Jidx)

class HeisenbergMeanFieldCalculator:
  def __init__(self):
    pass
  def set_site_resolved_magnetization(self, mag):
    self.sitemag = np.zeros((self.nsites), dtype=float) + np.array(mag, dtype=float)[:self.nsites]
  def set_magnetic_field(self, field):
    self.magfield = float(field)
  def set_exchange_parameters(self, param):
    self.exchange_parameters = np.array(param, dtype=float)
  def read_modelfile(self,infilename):
    infilehandle = open(infilename, 'r')
    lines = infilehandle.readlines()[2:]
    infilehandle.close()
    self.HeisenbergBonds = []
    self.MeanFieldBonds = []
    for i,l in enumerate(lines):
      sl = l.strip().split()
      if('#' not in l and len(sl) > 0):
        b = HeisenbergBond(sl[0], sl[1], sl[2])
        self.HeisenbergBonds.append(b)
      else:
        break
    for l in lines[i+3:]:
      sl = l.strip().split()
      if(len(sl) > 0):
        b = MeanFieldBond(sl[0], sl[1], sl[2])
        self.MeanFieldBonds.append(b)
    #determine number of sites from bonds
    self.nsites = max([b.s1 for b in self.HeisenbergBonds] + [b.s2 for b in self.HeisenbergBonds])+1
  def read_spinbasis(self, infilename):
    self.spinbasisfilename = infilename
    self.basisstates = []
    infilehandle = open(self.spinbasisfilename, 'r')
    for i,l in enumerate(infilehandle):
      if i>0:
        sl = l.strip().split()
        s = SpinState(sl[1])
        self.basisstates.append(s)
    infilehandle.close()
    self.nbasisstates = len(self.basisstates)
  def read_hamiltonian(self, infilename):
    self.hamiltonianmatrices = []
    infilehandle1 = open(infilename, 'r')
    for f in infilehandle1:
      infilehandle2 = open(f.strip(), 'r')
      i = []
      j = []
      v = []
      for l in infilehandle2:
        sl = l.strip().split()
        iidx = int(sl[0])
        jidx = int(sl[1])
        vval = float(sl[2])
        i.append(iidx)
        j.append(jidx)
        v.append(vval)
        #also add the lower triangle
        if iidx != jidx:
          i.append(jidx)
          j.append(iidx)
          v.append(vval)
      infilehandle2.close()
      #construct sparse matrix in coo format and convert directly to csr format
      self.hamiltonianmatrices.append(spsparse.coo_matrix((v,(i,j)), shape=(self.nbasisstates, self.nbasisstates)).tocsr())
    infilehandle1.close()
  def solve_cluster(self):
    #construct mean-field magnetic fields on each site
    sitefields = np.zeros((self.nsites), dtype=float)
    for b in self.MeanFieldBonds:
      sitefields[b.sr] -= self.exchange_parameters[b.Jidx]*self.sitemag[b.sm]
    #construct hamiltonian first as empty matrix
    A = spsparse.coo_matrix((self.nbasisstates, self.nbasisstates)).tocsr()
    #add hamiltonian matrices multiplied by respective exchange term
    for J,B in zip(self.exchange_parameters, self.hamiltonianmatrices):
      A = A + B.multiply(J)
    #add hamiltonian part from effective magnetic fields, this is a diagonal matrix by definition
    magterms = []
    i = []
    for j,s in enumerate(self.basisstates):
      i.append(j) 
      magterms.append(np.sum(np.multiply(s.get_sz_site_resolved(), sitefields)))
    M = spsparse.coo_matrix((magterms,(i,i)), shape=(self.nbasisstates, self.nbasisstates)).tocsr()
    A = A - M
    #initialize eigensolver  
    k = 1 #find only ground state
    self.eigenvalues, self.eigenvectors = splinalg.eigsh(A, k=k, which='SA')
  def calculate_new_magnetization(self):
    self.newsitemag = np.zeros(self.nsites)
    for s,p in zip(self.basisstates, self.eigenvectors):
      self.newsitemag += pow(p,2)*s.get_sz_site_resolved()
  def calculate_new_energy_per_site(self):
    mfenergy = 0.0
    for b in self.MeanFieldBonds:
      mfenergy -= 0.5*self.exchange_parameters[b.Jidx]*self.newsitemag[b.sr]*self.newsitemag[b.sm]
    self.energy_per_site = (self.eigenvalues[0] + mfenergy)/float(self.nsites)
  def get_magnetization_difference_measure_and_mix_magnetization(self):
    convparam = np.linalg.norm(self.newsitemag - self.sitemag)/self.nsites
    mixingfac = 0.45
    self.sitemag = self.sitemag*(1-mixingfac) + self.newsitemag*mixingfac
    return convparam
  def solve_selfconsistently(self):
    threshold = 1e-5
    maxiter = 1000
    itcounter = 0
    while True:
      self.solve_cluster()
      self.calculate_new_magnetization()
      self.calculate_new_energy_per_site()
      #print "Sitemag:", self.sitemag
      convparam = self.get_magnetization_difference_measure_and_mix_magnetization()
      #print 'Iteration: %i\nConvergence parameter: %f' % (itcounter, convparam)
      #print 'GS energy: % f' % self.energy_per_site
      itcounter += 1
      if(convparam < threshold):
        print 'Self-consistent loop converged after %i iterations.' % itcounter
        break
      if(itcounter == maxiter):
        print 'Self-consistent loop did not converge.'
        break
  def get_energy_per_site(self):
    return self.energy_per_site
  def get_total_magnetization_per_site(self):
    return np.mean(self.sitemag)
  def get_magnetization_of_site(self,idx):
    return self.sitemag[int(idx)]
  
class ClassicalGroundStateCalculator:
  def __init__(self):
    pass
  def read_modelfile(self,infilename):
    infilehandle = open(infilename, 'r')
    lines = infilehandle.readlines()[2:]
    infilehandle.close()
    self.HeisenbergBonds = []
    self.MeanFieldBonds = []
    for i,l in enumerate(lines):
      sl = l.strip().split()
      if('#' not in l and len(sl) > 0):
        b = HeisenbergBond(sl[0], sl[1], sl[2])
        self.HeisenbergBonds.append(b)
      else:
        break
    for l in lines[i+3:]:
      sl = l.strip().split()
      if(len(sl) > 0):
        b = MeanFieldBond(sl[0], sl[1], sl[2])
        self.MeanFieldBonds.append(b)
    #determine number of sites from bonds
    self.nsites = max([b.s1 for b in self.HeisenbergBonds] + [b.s2 for b in self.HeisenbergBonds])+1
  def set_exchange_parameters(self, param):
    self.exchange_parameters = np.array(param, dtype=float)
  def unique(self, iterable):
    seen = set()
    for x in iterable:
        if x in seen:
            continue
        seen.add(x)
        yield x
  def calculate_classical_energies(self):
    #first generate all unique possible states
    baseconfig = 5*[1] + 4*[-1]
    self.states = [s for s in self.unique(itertools.permutations(baseconfig))]
    self.energies = [self.get_energy_of_config_per_site(s) for s in self.states]
    self.sortidx = np.argsort(self.energies)
  def get_energy_of_config_per_site(self, state):
    energy = 0
    for b in self.HeisenbergBonds:
      energy += 0.25*state[b.s1]*state[b.s2]*self.exchange_parameters[b.Jidx]
    for b in self.MeanFieldBonds:
      energy += 0.125*state[b.sr]*state[b.sm]*self.exchange_parameters[b.Jidx] #half of normal bond, because it belongs to two unit cells
    return energy/float(self.nsites)
  def print_states_and_energies_per_site(self):
    for idx in self.sortidx[::-1]:
      print self.energies[idx], self.states[idx]
  def print_lowest_energy_and_state(self):
    idx = self.sortidx[0]
    print self.energies[idx], self.states[idx]