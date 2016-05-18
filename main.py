import os,sys 
from HeisenbergMeanField import HeisenbergMeanFieldCalculator

def main():
  modelfilename = 'examples/meanfieldmodel.dat'
  calc = HeisenbergMeanFieldCalculator()
  calc.set_magnetic_field(1.0)
  calc.read_modelfile(modelfilename)
  calc.solve_selfconsistently()
  return 0
  
main()