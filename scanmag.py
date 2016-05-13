import os
import subprocess
import numpy as np

def main():
  fieldvalues = np.linspace(0,1.5,num=20)
  chainlengths = range(2,13,2)
  fulldiag = 0
  periodic = 1
  J = -1.0
  nev = 1
  tempfilename = "test.dat"
  
  for L in chainlengths:
    magvals = []
    for B in fieldvalues:
      command = "./heisenberg %i %i %f %f %i %i" % (L, periodic, J, B, nev, fulldiag)
      p = subprocess.Popen(command.split(), stdout=subprocess.PIPE) 
      out, err = p.communicate()
      p.wait() 
      M = float(out[out.find("Bohr")+5:])
      magvals.append(M)
      
    outfilehandle = open("mag_L=%i_J=%1.2f.dat" % (L, J), 'w')
    outfilehandle.write("#B in T, <M> per site in Bohr\n")
    for B, M in zip(fieldvalues, magvals):
      outfilehandle.write("%f %f\n" % (B, M))
    outfilehandle.close()
      
 
main()

