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

struct HeisenbergChainParameters{
  int nsites, periodic, neigenvalues;
  fptype J, B, g;
};

class SpinState{
  public:
    SpinState();
    ~SpinState();
    SpinState* apply_splus(const int i);
    SpinState* apply_sminus(const int i);
    SpinState* apply_sz(const int i);
    fptype dot(SpinState* otherstate);
    void copy_to(SpinState* newstate);
    vector<bool> statevec;
    fptype prefactor, spin;
    bool copiedstate;
  private:

};

class SpinBasisGenerator{
  public:
    SpinBasisGenerator();
    ~SpinBasisGenerator();
    vector<SpinState> get_basis(const int nsites, const int sz);
  private:
    vector<bool> convert_to_binary_vectorbool(const int nsites, const unsigned long long int value);
    int get_sz(vector<bool> v);
};

class HeisenbergHamiltonianSolver{
  public:
    HeisenbergHamiltonianSolver();
    void set_bonds(vector<HeisenbergBond>& bonds);
    void set_fields(vector<SiteDependentMagneticField>& fields){this->fields = fields;};
    void calculate_eigenvalues_eigenvectors();
    void calculate_groundstate_site_dependent_magnetization();
  private:
    Eigen::MatrixXd get_hamiltonian(vector<SpinState>& basis);
    fptype get_hamiltonian_element(SpinState* u, SpinState* v);
    uint nsites;
    vector<HeisenbergBond> bonds;
    vector<SiteDependentMagneticField> fields;
    vector<int> allowed_sz;
    vector<vector<SpinState> > basis_sectors; 
    vector<Eigen::VectorXd> evals;
    vector<Eigen::MatrixXd> evecs;
    vector<fptype> maggs; //ground state site dependent magnetization
};

#endif

vector<int> get_allowed_sz(const int nsites);