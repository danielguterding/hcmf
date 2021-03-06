import os,sys 
import numpy as np
from HeisenbergMeanField import HeisenbergMeanFieldCalculator

class MagneticOrder:
  def __init__(self, name, pattern):
    self.name = str(name)
    self.pattern = np.array(pattern, dtype=float)

class MeanFieldResult:
  def __init__(self, J2, ordername, energypersite, totalmag, sitemag):
    self.index = float(J2)
    self.ordername = str(ordername)
    self.energypersite = float(energypersite)
    self.totalmag = float(totalmag)
    self.sitemag = float(sitemag)

def get_patterns():
  patterns = []
  patterns.append(MagneticOrder('Neel', [0.0, 0.5, -0.5, -0.5, 0.5, 0.0, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, 0.5, -0.5, 0.0, -0.5, 0.5, 0.0]))
  patterns.append(MagneticOrder('PM', 18*[0]))
  return patterns

def write_results(outfilename, results):
  outfilehandle = open(outfilename, 'w')
  headerstr = '#J2/J1'
  for entry in results[0]:
    headerstr += ', %s order energy/site, %s order moment/site ' % (entry.ordername, entry.ordername)
  headerstr += '\n'
  outfilehandle.write(headerstr)
  
  for res_set in results:
    outstr = '% 1.7f ' % res_set[0].index 
    for entry in res_set:
      outstr += '% 1.7f %1.7f ' % (entry.energypersite, abs(entry.sitemag))
    outstr += '\n'
    outfilehandle.write(outstr)
  outfilehandle.close()

def main():
  outfilename = 'results_Y3Cu9OH19Cl8_18spin_cluster2_J0_0p1.dat'
  modelfilename = 'examples/meanfieldmodel_Y3Cu9OH19Cl8_18spin_cluster2.dat'
  spinbasisfilename = 'spinbasis_Y3Cu9OH19Cl8_18spin_cluster2.dat'
  hamiltonianfilename = 'hamiltonian_Y3Cu9OH19Cl8_18spin_cluster2.dat'
  magneticpatterns = get_patterns()
  calc = HeisenbergMeanFieldCalculator()
  calc.read_modelfile(modelfilename)
  calc.read_spinbasis(spinbasisfilename)
  calc.read_hamiltonian(hamiltonianfilename)
  J0 = 0.1
  J1 = 1.0
  J2vals = np.linspace(0.1,1.0,num=10,endpoint=True)
  B = 0.0
  results = []
  for J2 in J2vals:
    res_this = []
    for p in magneticpatterns:
      exchange_parameters = [J0, J1, J2]
      calc.set_magnetic_field(B)
      calc.set_exchange_parameters(exchange_parameters)
      calc.set_site_resolved_magnetization(p.pattern)
      calc.solve_selfconsistently()
      res_this.append(MeanFieldResult(J2, p.name,calc.get_energy_per_site(),calc.get_total_magnetization_per_site(),calc.get_magnetization_of_site(0)))
    results.append(res_this)
  write_results(outfilename, results)
  return 0
  
main()