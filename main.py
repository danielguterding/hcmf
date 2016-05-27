import os,sys 
import numpy as np
from HeisenbergMeanField import HeisenbergMeanFieldCalculator

def main():
  modelfilename = 'examples/meanfieldmodel.dat'
  exchange_parameters = [1.0, 0.3]
  calc = HeisenbergMeanFieldCalculator()
  calc.read_modelfile(modelfilename)
  calc.set_site_resolved_magnetization([-1, 1, 1, -1, -1, 1, 1, -1])
  #for B in np.linspace(0,1,num=100):
  for B in [0]:
    calc.set_magnetic_field(B)
    calc.set_exchange_parameters(exchange_parameters)
    calc.solve_selfconsistently()
    print B,calc.get_energy_per_site(),calc.get_total_magnetization_per_site()
  return 0
  
main()