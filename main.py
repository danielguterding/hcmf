import os,sys 
import numpy as np
from HeisenbergMeanField import HeisenbergMeanFieldCalculator

def main():
  modelfilename = 'examples/meanfieldmodel.dat'
  calc = HeisenbergMeanFieldCalculator()
  calc.read_modelfile(modelfilename)
  #for B in np.linspace(0,1,num=100):
  for B in [0]:
    calc.set_magnetic_field(B)
    calc.solve_selfconsistently()
    print B,calc.get_energy_per_site(),calc.get_total_magnetization_per_site()
  return 0
  
main()