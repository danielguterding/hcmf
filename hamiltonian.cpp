//hamiltonian.cpp

#include "hamiltonian.hpp"

SpinBasisGeneratorSZ::SpinBasisGeneratorSZ(const uint nsites, const int sz){
  //sz is the net spin is measured in units of 1/2 to make basis generation easier
  //this->basis.resize(0);
  
  const uint nup = nsites/2 + sz;
  //determine minimum integer value correspoding to bitpattern with all up spins to left and right
  unsigned long long int maxbin = 0;
   unsigned long long int minbin = 0;
  for(uint i=0;i<nup;i++){
    maxbin += pow(2, nsites-i-1); 
    minbin += pow(2, i);
  }
  maxbin++;
  
  string bitstring;
  for(unsigned long long int i=minbin;i<maxbin;i++){
    vector<bool> b = convert_to_binary_vectorbool(nsites, i);
    if(sz == get_sz(b)){
      SpinState s;
      s.statevec = b;
      basis.push_back(s);
    }
  }
}

SpinBasisGeneratorSZ::~SpinBasisGeneratorSZ(){
  
}

vector<bool> SpinBasisGeneratorSZ::convert_to_binary_vectorbool(const int nsites, const unsigned long long int value){
  
  vector<bool> vec(nsites);
  fill(vec.begin(), vec.end(), 0);

  for (int current_bit = nsites-1; current_bit >= 0; current_bit--){ //does not support more than 64 sites
    if ((value & (1ULL << current_bit)) != 0){
      vec[current_bit] = 1;        
    }
    else{
      vec[current_bit] = 0;      
    }
  }

  return vec;
}

int SpinBasisGeneratorSZ::get_sz(vector<bool> v){
  
  int sz = 0;
  for(uint i=0;i<v.size();i++){
    sz += 2*v[i]-1; //substract one if entry is zero, add one if entry is one
  }
  
  return sz;
}

HeisenbergHamiltonianCalculator::HeisenbergHamiltonianCalculator(){
  
}

void HeisenbergHamiltonianCalculator::set_bonds(vector<HeisenbergBond>& bonds){
  
  this->bonds = bonds;
  sort(this->bonds.begin(), this->bonds.end(), [](const HeisenbergBond& b1, const HeisenbergBond& b2){return b1.s2 > b2.s2;}); //compare with inline lambda function
  this->nsites = this->bonds[0].s2+1;
  
  //sort bonds by Jidx to determine how many different interaction terms we have
  sort(this->bonds.begin(), this->bonds.end(), [](const HeisenbergBond& b1, const HeisenbergBond& b2){return b1.Jidx > b2.Jidx;}); //compare with inline lambda function
  this->nJidx = bonds[0].Jidx+1;
}

void HeisenbergHamiltonianCalculator::set_basis(const vector<SpinState>& basis){
  
  this->basis = basis;
}

void HeisenbergHamiltonianCalculator::set_outfilename(const string outfilename){
  
  this->outfilename = outfilename;
  boost::filesystem::path outfilepath(this->outfilename);
  boost::filesystem::ofstream outfilehandle(outfilepath);
  this->outfilenameseachinteraction.resize(0);
  for(uint i=0;i<this->nJidx+1;i++){
    string ofn = (boost::format("%s.%03i") % this->outfilename % i).str();
    outfilenameseachinteraction.push_back(ofn);
    outfilehandle << ofn << endl;
  }
  outfilehandle.close();
}

void HeisenbergHamiltonianCalculator::calculate_elements(){
  
  const fptype threshold = 1e-12;
  HamiltonianElement element;
  for(uint i=0;i<bonds.size();i++){//loop over bonds
    boost::filesystem::path outfilepath(this->outfilenameseachinteraction[bonds[i].Jidx]);
    boost::filesystem::ofstream outfilehandle(outfilepath);
    for(unsigned long long int j=0;j<basis.size();j++){//row index
      for(unsigned long long int k=j;k<basis.size();k++){//column index
        element.i = j;
        element.j = k;
        element.v = get_hamiltonian_element(&bonds[i], &basis[j], &basis[k]);
        if(fabs(element.v) > threshold){
          outfilehandle << boost::format("%14i %14i % 1.14f") % element.i % element.j % element.v << endl;
        }
      }
    }
    outfilehandle.close();
  }
}

fptype HeisenbergHamiltonianCalculator::get_hamiltonian_element(const HeisenbergBond* b, SpinState* u, SpinState* v){
  
  fptype element = 0;
  element += 0.5*(v->dot(u->apply_sminus(b->s1)->apply_splus(b->s2)) + v->dot(u->apply_sminus(b->s2)->apply_splus(b->s1)));
  element += v->dot(u->apply_sz(b->s2)->apply_sz(b->s1));
  return element;
};