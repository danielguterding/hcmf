//hamiltonian.cpp

#include "hamiltonian.hpp"

SpinState::SpinState(){
  
  this->spin = 0.5;
  this->prefactor = 1.0;
  this->copiedstate = false;
}

SpinState::~SpinState(){
  
}

SpinState* SpinState::apply_splus(const int i){
  
  if(!this->copiedstate){
    SpinState* newstate = new SpinState;
    this->copy_to(newstate);
    if(!newstate->statevec[i]){ //if bit is set to 0
      newstate->statevec[i] = 1;
    }
    else{
      newstate->prefactor = 0;
    }
    return newstate;
  }
  else{
    if(!this->statevec[i]){ //if bit is set to 0
      this->statevec[i] = 1;
    }
    else{
      this->prefactor = 0;
    }
    return this;
  }
}

SpinState* SpinState::apply_sminus(const int i){
  
  if(!this->copiedstate){
    SpinState* newstate = new SpinState;
    this->copy_to(newstate);
    if(newstate->statevec[i]){ //if bit is set to 1
      newstate->statevec[i] = 0;
    }
    else{
      newstate->prefactor = 0;
    }
    return newstate;
  }
  else{
    if(this->statevec[i]){ //if bit is set to 1
      this->statevec[i] = 0;
    }
    else{
      this->prefactor = 0;
    }
    return this;
  }
}

SpinState* SpinState::apply_sz(const int i){
  
  if(!this->copiedstate){
    SpinState* newstate = new SpinState;
    this->copy_to(newstate);
    newstate->prefactor *= newstate->spin * (2*newstate->statevec[i]-1);
    return newstate;
  }
  else{
    this->prefactor *= this->spin * (2*this->statevec[i]-1);
    return this;
  }
}

fptype SpinState::dot(SpinState* otherstate){
  
  bool allequal = true;
  const uint nentries = this->statevec.size();
  uint i=0;
  while((i < nentries) && allequal){
    if(this->statevec[i] != otherstate->statevec[i]){
      allequal = false;
    }
    i++;
  }
  
  if(allequal){
    return this->prefactor*otherstate->prefactor;
  }
  else{
    return 0.0;
  }
}

void SpinState::copy_to(SpinState* newstate){
  
  newstate->statevec = this->statevec;
  newstate->prefactor = this->prefactor;
  newstate->spin = this->spin;
  newstate->copiedstate = true;
}

SpinBasisGenerator::SpinBasisGenerator(){
  
}

SpinBasisGenerator::~SpinBasisGenerator(){
  
}

vector<SpinState> SpinBasisGenerator::get_basis(const int nsites, const int sz){
  //sz is the net spin is measured in units of 1/2 to make basis generation easier
  
  vector<SpinState> basis;
  string bitstring;
  for(uint i=0;i<uint(pow(2,nsites));i++){
    vector<bool> b = convert_to_binary_vectorbool(nsites, i);
    if(sz == get_sz(b)){
      SpinState s;
      s.statevec = b;
      basis.push_back(s);
    }
  }
  
  return basis;
}

vector<bool> SpinBasisGenerator::convert_to_binary_vectorbool(const int nsites, const unsigned long long int value){
  
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

int SpinBasisGenerator::get_sz(vector<bool> v){
  
  int sz = 0;
  for(uint i=0;i<v.size();i++){
    sz += 2*v[i]-1; //substract one if entry is zero, add one if entry is one
  }
  
  return sz;
}

HeisenbergHamiltonianSolver::HeisenbergHamiltonianSolver(){
  
}

void HeisenbergHamiltonianSolver::set_bonds(vector<HeisenbergBond>& bonds){
  
  this->bonds = bonds;
  sort(this->bonds.begin(), this->bonds.end(), [](const HeisenbergBond& b1, const HeisenbergBond& b2){return b1.s2 < b2.s2;}); //compare with inline lambda function
  this->nsites = this->bonds[0].s2;
}

void HeisenbergHamiltonianSolver::calculate_eigenvalues_eigenvectors(){
  
  allowed_sz = get_allowed_sz(this->nsites); //get allowed sz sectors
  evals.resize(0);
  evecs.resize(0);
  basis_sectors.resize(0);
  SpinBasisGenerator bgen;
  for(uint i=0;i<allowed_sz.size();i++){
    basis_sectors.push_back(bgen.get_basis(this->nsites,allowed_sz[i]));

    Eigen::MatrixXd h = get_hamiltonian(basis_sectors[i]);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    solver.compute(h);
    evals.push_back(solver.eigenvalues());
    evecs.push_back(solver.eigenvectors());
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
  
  const fptype gfac = 2.0;
  fptype element = 0;
  for(uint i=0;i<bonds.size();i++){ //loop over bonds
    const HeisenbergBond b = bonds[i];
    element += -b.J/2.0*(v->dot(u->apply_sminus(b.s1)->apply_splus(b.s2)) + v->dot(u->apply_sminus(b.s2)->apply_splus(b.s1))) - b.J*v->dot(u->apply_sz(b.s2)->apply_sz(b.s1));
  }
  
  for(uint i=0;i<fields.size();i++){ //loop over magnetic field entries
    const SiteDependentMagneticField f = fields[i];
    element -= gfac*f.h*v->dot(u->apply_sz(f.s));
  }
  
  return element;
};

void HeisenbergHamiltonianSolver::calculate_groundstate_site_dependent_magnetization(){
  
  maggs.resize(0);
  //first identify ground state energy
  fptype gsenergy = evals[0](0);
  for(uint i=1;i<evals.size();i++){
    if(evals[i](0) < gsenergy){
      gsenergy = evals[i](0);
    }
  }
  //use states that lie in very narrow region around ground state and calculate their magnetization
  
}

vector<int> get_allowed_sz(const int nsites){
  //sz is the net spin measured in units of 1/2
  vector<int> sz;
  for(int i=-nsites;i<(nsites+1);i+=2){
    sz.push_back(i);
  }
  
  return sz;
}