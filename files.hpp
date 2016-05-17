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
  fptype J; //exchange coupling
};

struct SiteDependentMagneticField{
  uint s; //site index
  fptype h; //value of magnetic field on that site
};

class HeisenbergBondReader{
  public:
    HeisenbergBondReader();
    void read_file(string infilename);
    vector<HeisenbergBond> get_bonds(){return bonds;};
  private:
    vector<HeisenbergBond> bonds;
};

class SiteDependentMagneticFieldReader{
  public:
    SiteDependentMagneticFieldReader();
    void read_file(string infilename);
    vector<SiteDependentMagneticField> get_fields(){return fields;};
  private:
    vector<SiteDependentMagneticField> fields;
};

#endif

string trim_all(const std::string &str);