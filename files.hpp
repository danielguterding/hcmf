//files.hpp
#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include "typedefs.hpp"

using namespace std;

#ifndef __FILES_H_INCLUDED__
#define __FILES_H_INCLUDED__

struct HeisenbergBond{
  uint s1,s2; //site indices
  fptype Jidx; //exchange index
};

struct HamiltonianElement{
  unsigned long long int i,j; //matrix indices
  fptype v; //value
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

class HeisenbergBondReader{
  public:
    HeisenbergBondReader();
    void read_file(const string infilename);
    vector<HeisenbergBond> get_bonds(){return bonds;};
    uint get_nsites(){return nsites;};
  private:
    vector<HeisenbergBond> bonds;
    uint nsites;
};

class SpinBasisWriter{
  public:
    SpinBasisWriter();
    void write_basis(string outfilename, vector<SpinState> basis);
  private:
  
};

#endif

string trim_all(const std::string &str);