import os,sys 
from HeisenbergMeanField import HeisenbergMeanFieldCalculator

def main():
  modelfilename = 'examples/meanfieldmodel.dat'
  calc = HeisenbergMeanFieldCalculator()
  calc.read_modelfile(modelfilename)
  calc.solve_selfconsistently()
  return 0
  
main()