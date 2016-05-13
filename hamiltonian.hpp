//hamiltonian.hpp
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

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

class HeisenbergChainEvaluator{
  public:
    HeisenbergChainEvaluator(HeisenbergChainParameters& p);
    ~HeisenbergChainEvaluator();
    Eigen::MatrixXd get_hamiltonian_dense(vector<SpinState>& basis);
    Eigen::SparseMatrix<fptype> get_hamiltonian_sparse(vector<SpinState>& basis, fptype element_threshold=1e-7);
  private:
    fptype get_element(SpinState* u, SpinState* v);
    HeisenbergChainParameters p;
};

class LanczosSolver{
  public:
    LanczosSolver();
    ~LanczosSolver();
    void compute(Eigen::SparseMatrix<fptype> m, int nev, const fptype ev_threshold=1e-5);
    Eigen::VectorXd get_eigenvalues(){return eigenvalues;};
  private:
    Eigen::MatrixXd hessenberg;
    Eigen::VectorXd eigenvalues;
};

#endif

vector<int> get_allowed_sz(const int nsites);