//hamiltonian.hpp
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "files.hpp"
#include "typedefs.hpp"

using namespace std;

#ifndef __HAMILTONIAN_H_INCLUDED__
#define __HAMILTONIAN_H_INCLUDED__

class SpinBasisGeneratorSZ{
  public:
    SpinBasisGeneratorSZ(const uint nsites, const int sz);
    ~SpinBasisGeneratorSZ();
    vector<SpinState> get_basis(){return basis;};
  private:
    vector<bool> convert_to_binary_vectorbool(const int nsites, const unsigned long long int value);
    int get_sz(vector<bool> v);
    vector<SpinState> basis;
};

/*
class HeisenbergHamiltonianSolver{
  public:
    HeisenbergHamiltonianSolver();
    void set_bonds(vector<HeisenbergBond>& bonds);
    void calculate_eigenvalues_eigenvectors();
  private:
    Eigen::MatrixXd get_hamiltonian(vector<SpinState>& basis);
    fptype get_hamiltonian_element(SpinState* u, SpinState* v);
    uint nsites;
    fptype gsenergy, totalmag;
    vector<HeisenbergBond> bonds;
    vector<int> allowed_sz;
    vector<vector<SpinState> > basis_sectors; 
};
*/

#endif

vector<int> get_allowed_sz(const int nsites);