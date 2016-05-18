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
    const string InfilenameHeisenbergBonds = argv[1];
    const string InfilenameSiteDependentMagneticFields = argv[2];
    const string OutfilenameMagnetization = argv[3];
    
    //read Heisenberg bonds
    HeisenbergBondReader BondReader;
    BondReader.read_file(InfilenameHeisenbergBonds);
    vector<HeisenbergBond> bonds = BondReader.get_bonds();
    
    //read magnetic field on each site
    SiteDependentMagneticFieldReader FieldReader;
    FieldReader.read_file(InfilenameSiteDependentMagneticFields);
    vector<SiteDependentMagneticField> fields = FieldReader.get_fields();
    
    //solve Hamiltonian
    HeisenbergHamiltonianSolver solver;
    solver.set_bonds(bonds);
    solver.set_fields(fields);
    solver.calculate_eigenvalues_eigenvectors();
    solver.calculate_groundstate_site_dependent_magnetization();
    
    //get magnetization per site and write to file
    const fptype totalmag = solver.get_groundstate_total_magnetization_per_site();
    const vector<fptype> maggs = solver.get_groundstate_site_dependent_magnetization();
    const fptype energy_per_site = solver.get_groundstate_energy_per_site();
    MagnetizationWriter MagWriter;
    MagWriter.set_energy_per_site(energy_per_site);
    MagWriter.set_total_magnetization(totalmag);
    MagWriter.set_site_resolved_magnetization(maggs);
    MagWriter.write_magnetization(OutfilenameMagnetization);
    cout << "Output written. Program finished successfully." << endl;
    return 0;
  }
  else{
    cout << "Wrong number of input arguments. Please supply input file names for Heisenberg Bonds and site-dependent magnetic fields and output file name for site-dependent magnetization." << endl;
    return 1;
  }
}