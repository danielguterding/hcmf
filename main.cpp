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
  
  if(3 == argc){
    const string InfilenameBonds = argv[1];
    const string OutfilenameSpinBasis = argv[2];
    
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

    //solve Hamiltonian
    /* HeisenbergHamiltonianSolver solver;
    solver.set_bonds(bonds);
    solver.calculate_eigenvalues_eigenvectors(); */
    
    cout << "Output written. Program finished successfully." << endl;
    return 0;
  }
  else{
    cout << "Wrong number of input arguments. Please supply input file names for bonds." << endl;
    return 1;
  }
}