import os
import subprocess
import numpy as np

class HeisenbergBond:
  def __init__(self, s1, s2, J):
    self.s1 = int(s1)
    self.s2 = int(s2)
    self.J = float(J)

class MeanFieldBond:
  def __init__(self, s1, s2, J):
    self.sr = int(s1) #real site 
    self.sm = int(s2) #mirror site outside of cluster
    self.J = float(J)

class HeisenbergMeanFieldCalculator:
  def __init__(self):
    self.filename_sitemag = 'sitemag.dat'
    self.filename_fields = 'fields.dat'
    self.filename_heisenberg_bonds = 'bonds.dat'
    self.magfield = 0.0
  def set_magnetic_field(self, field):
    self.magfield = float(field)
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
    self.sitemag = np.zeros((self.nsites), dtype=float)
  def write_input_for_cluster_solver(self):
    #write Heisenberg bonds input file for cluster solver
    outfilehandle = open(self.filename_heisenberg_bonds, 'w')
    outfilehandle.write('#idx site 1, idx site 2, exchange coupling J\n')
    for b in self.HeisenbergBonds:
      outfilehandle.write('%i %i % f\n' % (b.s1, b.s2, b.J))
    outfilehandle.close()
    #write effective magnetic fields for cluster solver
    outfilehandle = open(self.filename_fields, 'w')
    outfilehandle.write('#site idx, magnetic field\n')
    sitefields = self.magfield*np.ones((self.nsites), dtype=float)
    for b in self.MeanFieldBonds: #add mean field bond induced local fields to global magnetic field
      sitefields[b.sr] -= b.J*self.sitemag[b.sm]
    for i,f in enumerate(sitefields):
      outfilehandle.write('%i % f\n' % (i, f))
    outfilehandle.close()
  def execute_cluster_solver(self):
    command = './heisenberg %s %s %s' % (self.filename_heisenberg_bonds, self.filename_fields, self.filename_sitemag)
    FNULL = open(os.devnull, 'w')
    p = subprocess.Popen(command.split(), stdout=FNULL, stderr=subprocess.STDOUT) #hide cluster solver output
    p.wait()
    FNULL.close()
  def read_site_magnetization(self):
    infilehandle = open(self.filename_sitemag)
    lines = infilehandle.readlines()
    #self.totalmag_per_site = float(lines[1])
    self.newsitemag = np.zeros((self.nsites))
    for l in lines[3:]:
      sl = l.strip().split()
      i = int(sl[0])
      h = float(sl[1])
      self.newsitemag[i] = h
    infilehandle.close()
  def get_magnetization_difference_measure_and_mix_magnetization(self):
    convparam = np.linalg.norm(self.newsitemag - self.sitemag)/self.nsites
    mixingfac = 0.4
    self.sitemag = self.sitemag*(1-mixingfac) + self.newsitemag*mixingfac
    return convparam
  def solve_selfconsistently(self):
    threshold = 1e-3
    itcounter = 0
    while True:
      self.write_input_for_cluster_solver()
      self.execute_cluster_solver()
      self.read_site_magnetization()
      convparam = self.get_magnetization_difference_measure_and_mix_magnetization()
      print 'Iteration: %i\nConvergence parameter: %f' % (itcounter, convparam)
      itcounter += 1
      if(convparam < threshold):
        print 'Self-consistent loop converged after %i iterations.' % itcounter
        break
        