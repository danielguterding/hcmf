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

def get_patterns_honeycomb():
  patterns = []
  patterns.append(MagneticOrder('PM', [0, 0, 0, 0, 0, 0, 0, 0, 0]))
  #patterns.append(MagneticOrder('Trimer', [0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5]))
  #patterns.append(MagneticOrder('Hexagon Neel', [-0.5, -0.5, -0.5, 0.5, 0.5, -0.5, 0.5, -0.5, 0.5]))
  return patterns

def write_results(outfilename, results):
  outfilehandle = open(outfilename, 'w')
  headerstr = '#J2/J3'
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
  outfilename = 'result_meanfieldmodel_Y3Cu9OH19Cl8_9spin_cluster3.dat'
  modelfilename = 'examples/meanfieldmodel_Y3Cu9OH19Cl8_9spin_cluster3.dat'
  magneticpatterns = get_patterns_honeycomb()
  calc = HeisenbergMeanFieldCalculator()
  calc.read_modelfile(modelfilename)
  J2vals = np.linspace(0.0, 1,num=40,endpoint=True)
  B = 0.0
  results = []
  for J2 in J2vals:
    res_this_J2 = []
    for p in magneticpatterns:
      exchange_parameters = [0.3, J2, 1.0]
      calc.set_magnetic_field(B)
      calc.set_exchange_parameters(exchange_parameters)
      calc.set_site_resolved_magnetization(p.pattern)
      calc.solve_selfconsistently()
      res_this_J2.append(MeanFieldResult(J2, p.name,calc.get_energy_per_site(),calc.get_total_magnetization_per_site(),calc.get_magnetization_of_site(0)))
    results.append(res_this_J2)
  write_results(outfilename, results)
  return 0
  
main()