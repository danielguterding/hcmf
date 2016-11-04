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

/*

HeisenbergHamiltonianSolver::HeisenbergHamiltonianSolver(){
  
}

void HeisenbergHamiltonianSolver::set_bonds(vector<HeisenbergBond>& bonds){
  
  this->bonds = bonds;
  sort(this->bonds.begin(), this->bonds.end(), [](const HeisenbergBond& b1, const HeisenbergBond& b2){return b1.s2 > b2.s2;}); //compare with inline lambda function
  this->nsites = this->bonds[0].s2+1;
}

void HeisenbergHamiltonianSolver::calculate_eigenvalues_eigenvectors(){
  
  allowed_sz = get_allowed_sz(this->nsites); //get allowed sz sectors
  evals.resize(0);
  evecs.resize(0);
  basis_sectors.resize(0);
  SpinBasisGenerator bgen;
  //solve Hamiltonian for each sector of sz
  for(uint i=0;i<allowed_sz.size();i++){
    basis_sectors.push_back(bgen.get_basis(this->nsites,allowed_sz[i]));

    Eigen::MatrixXd h = get_hamiltonian(basis_sectors[i]);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    solver.compute(h);
    evals.push_back(solver.eigenvalues());
    //cout << evals[i].transpose() << endl;
    evecs.push_back(solver.eigenvectors());
  }
  //identify ground state energy
  gsenergy = evals[0](0);
  for(uint i=1;i<evals.size();i++){
    if(evals[i](0) < gsenergy){
      gsenergy = evals[i](0);
    }
  }
}

Eigen::MatrixXd HeisenbergHamiltonianSolver::get_hamiltonian(vector<SpinState>& basis){
  
  const int nstates = basis.size();
  Eigen::MatrixXd h = Eigen::MatrixXd::Zero(nstates, nstates);
  fptype element=0;
  for(int i=0;i<nstates;i++){
    for(int j=i;j<nstates;j++){ //only calculate upper half of the matrix, because it is real symmetric
      element = this->get_hamiltonian_element(&basis[i], &basis[j]);
      if(i == j){
        h(i,i) = element;
      }
      else{
        h(i,j) = element;
        h(j,i) = h(i,j);
      }
    }
  }
  
  return h;
}

fptype HeisenbergHamiltonianSolver::get_hamiltonian_element(SpinState* u, SpinState* v){
  
  fptype element = 0;
  for(uint i=0;i<bonds.size();i++){ //loop over bonds
    const HeisenbergBond b = bonds[i];
    element += b.J/2.0*(v->dot(u->apply_sminus(b.s1)->apply_splus(b.s2)) + v->dot(u->apply_sminus(b.s2)->apply_splus(b.s1)));
    element += b.J*v->dot(u->apply_sz(b.s2)->apply_sz(b.s1));
  }
  
  return element;
};

*/

vector<int> get_allowed_sz(const int nsites){
  //sz is the net spin measured in units of 1/2
  vector<int> sz;
  for(int i=-nsites;i<(nsites+1);i+=2){
    sz.push_back(i);
  }
  
  return sz;
}