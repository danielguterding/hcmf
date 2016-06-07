import os
import subprocess
import numpy as np
import itertools

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
    self.filename_sitemag = 'sitemag.dat'
    self.filename_fields = 'fields.dat'
    self.filename_heisenberg_bonds = 'bonds.dat'
    self.magfield = 0.0
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
  def write_input_for_cluster_solver(self):
    #write Heisenberg bonds input file for cluster solver
    outfilehandle = open(self.filename_heisenberg_bonds, 'w')
    outfilehandle.write('#idx site 1, idx site 2, exchange coupling J\n')
    for b in self.HeisenbergBonds:
      outfilehandle.write('%i %i % f\n' % (b.s1, b.s2, self.exchange_parameters[b.Jidx]))
    outfilehandle.close()
    #write effective magnetic fields for cluster solver
    outfilehandle = open(self.filename_fields, 'w')
    outfilehandle.write('#site idx, magnetic field\n')
    sitefields = 2.0*self.magfield*np.ones((self.nsites), dtype=float)
    for b in self.MeanFieldBonds: #add mean field bond induced local fields to global magnetic field
      addterm = -self.exchange_parameters[b.Jidx]*self.sitemag[b.sm]
      #print addterm
      sitefields[b.sr] += addterm 
    #print sitefields
    for i,f in enumerate(sitefields):
      outfilehandle.write('%i % f\n' % (i, f))
    outfilehandle.close()
  def execute_cluster_solver(self):
    command = './heisenberg %s %s %s' % (self.filename_heisenberg_bonds, self.filename_fields, self.filename_sitemag)
    FNULL = open(os.devnull, 'w')
    p = subprocess.Popen(command.split(), stdout=FNULL, stderr=subprocess.STDOUT) #hide cluster solver output
    #p = subprocess.Popen(command.split())
    p.wait()
    FNULL.close()
  def read_site_magnetization(self):
    infilehandle = open(self.filename_sitemag)
    lines = infilehandle.readlines()
    self.energy_per_site = float(lines[1])
    self.totalmag_per_site = float(lines[3])
    self.newsitemag = np.zeros((self.nsites))
    for l in lines[5:]:
      sl = l.strip().split()
      i = int(sl[0])
      h = float(sl[1])
      self.newsitemag[i] = h
    #add mean field terms <Si><Sj> to energy
    addenergy = 0.0
    for b in self.MeanFieldBonds:
      addenergy -= 0.5*self.exchange_parameters[b.Jidx]*self.newsitemag[b.sr]*self.newsitemag[b.sm]
    self.energy_per_site += addenergy/float(self.nsites)
    infilehandle.close()
  def get_magnetization_difference_measure_and_mix_magnetization(self):
    #print self.sitemag
    #print self.newsitemag
    convparam = np.linalg.norm(self.newsitemag - self.sitemag)/self.nsites
    mixingfac = 0.45
    self.sitemag = self.sitemag*(1-mixingfac) + self.newsitemag*mixingfac
    return convparam
  def solve_selfconsistently(self):
    threshold = 1e-5
    maxiter = 1000
    itcounter = 0
    while True:
      self.write_input_for_cluster_solver()
      self.execute_cluster_solver()
      self.read_site_magnetization()
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
    return self.totalmag_per_site
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
    #frist generate all unique possible states
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