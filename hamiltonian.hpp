//hamiltonian.hpp
#include <iostream>
#include <vector>

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

class HeisenbergHamiltonianCalculator{
  public:
    HeisenbergHamiltonianCalculator();
    void set_bonds(vector<HeisenbergBond>& bonds);
    void set_basis(const vector<SpinState>& basis);
    void set_outfilename(const string outfilename);
    void calculate_elements();
  private:
    fptype get_hamiltonian_element(const HeisenbergBond* b, SpinState* u, SpinState* v);
    uint nsites, nJidx;
    vector<HeisenbergBond> bonds;
    vector<SpinState> basis;
    string outfilename;
};

#endif

vector<int> get_allowed_sz(const int nsites);