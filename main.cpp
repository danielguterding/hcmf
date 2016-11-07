//main.cpp
#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "files.hpp"
#include "hamiltonian.hpp"
#include "typedefs.hpp"

using namespace std;

int main(int argc, char* argv[]){
  
  if(4 == argc){
    const string InfilenameBonds = argv[1];
    const string OutfilenameSpinBasis = argv[2];
    const string OutfilenameHamiltonian = argv[3];
    
    //read Heisenberg bonds
    HeisenbergBondReader HeisenbergBondReader;
    HeisenbergBondReader.read_file(InfilenameBonds);
    vector<HeisenbergBond> bonds = HeisenbergBondReader.get_bonds();

    //get number of sites
    const uint nsites = HeisenbergBondReader.get_nsites();
    //set desired dz value to zero
    const int sz = 0;
    
    //get and write spin basis
    SpinBasisGeneratorSZ BasisGenerator(nsites, sz);
    vector<SpinState> basis = BasisGenerator.get_basis();
    
    //write basis to disk
    SpinBasisWriter BasisWriter;
    BasisWriter.write_basis(OutfilenameSpinBasis, basis);
    
    //calculate Hamiltonian
    HeisenbergHamiltonianCalculator calc;
    calc.set_bonds(bonds);
    calc.set_outfilename(OutfilenameHamiltonian);
    calc.set_basis(basis);
    calc.calculate_elements();

    //solve Hamiltonian
    /* HeisenbergHamiltonianSolver solver;
    solver.set_bonds(bonds);
    solver.calculate_eigenvalues_eigenvectors(); */
    
    cout << "Output written. Program finished successfully." << endl;
    return 0;
  }
  else{
    cout << "Wrong number of input arguments. Please supply input file names for bonds, output file name for the spin basis and output file name for the Hamiltonian." << endl;
    return 1;
  }
}