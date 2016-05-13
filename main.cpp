//main.cpp
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "hamiltonian.hpp"
#include "typedefs.hpp"

using namespace std;

int main(int argc, char* argv[]){
  
  if(7 == argc){
    HeisenbergChainParameters p;
    p.nsites = atoi(argv[1]);
    p.periodic = atoi(argv[2]);
    p.J = atof(argv[3]);
    p.B = atof(argv[4]);
    p.g = 2.0;
    p.neigenvalues = atoi(argv[5]);
    const bool fulldiag = atoi(argv[6]);
    
    cout << boost::lexical_cast<string>(boost::format("Diagonalization program for the Heisenberg chain (L=%i, J=%1.2f).") % p.nsites % p.J) << endl;
    
    if(fulldiag){
      cout << "Using full diagonalization of the Hamiltonian." << endl;
    }
    else{
      cout << "Using Lanczos diagonalization of the Hamiltonian." << endl;
    }
    
    if(p.periodic){
      cout << "Periodic boundary conditions enabled." << endl;
    }
    else{
      cout << "Periodic boundary conditions disabled." << endl;
    }
    
    vector<int> allowed_sz = get_allowed_sz(p.nsites); //get allowed sz sectors
    SpinBasisGenerator bgen;
    HeisenbergChainEvaluator heisenbergchaineval(p);
    fptype groundstateenergy=0;
    int groundstatesz = 0;
    bool groundstateenergyset = false;
    for(uint i=0;i<allowed_sz.size();i++){
      cout << boost::lexical_cast<string>(boost::format("Solution for N=%i, Sz=%1.1f:") % p.nsites % allowed_sz[i]) << endl;
      
      vector<SpinState> basis = bgen.get_basis(p.nsites,allowed_sz[i]);
      fptype ev0 = 0;
      Eigen::VectorXd eigenvalues;
      
      if(fulldiag || (basis.size() < 3)){
	Eigen::MatrixXd h = heisenbergchaineval.get_hamiltonian_dense(basis);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
	solver.compute(h);
	eigenvalues = solver.eigenvalues();
	ev0 = eigenvalues(0);
      }
      else{
	Eigen::SparseMatrix<fptype> h = heisenbergchaineval.get_hamiltonian_sparse(basis);
	LanczosSolver solver;
	solver.compute(h, p.neigenvalues);
	eigenvalues = solver.get_eigenvalues();
	ev0 = eigenvalues(0);
      }
      
      if((!groundstateenergyset) || (ev0 < groundstateenergy)){
	groundstateenergy = ev0;
	groundstatesz = allowed_sz[i];
	groundstateenergyset = true;
      }
      
      cout << "Eigenvalues:" << endl << eigenvalues.head(min(p.neigenvalues, int(eigenvalues.rows()))).transpose() << endl << endl;
    }
    
    fptype magnetization_per_site = groundstatesz/float(p.nsites);
    cout << "Magnetization per site in Bohr: " << magnetization_per_site << endl;
  
  }
  else{
    cout << "Wrong number of input arguments. Please supply number of sites, periodic boundary switch, Heisenberg J, magnetic field, number of eigenvalues and full diagonalization switch." << endl;
  }
  
  return 0;
}