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

HeisenbergChainEvaluator::HeisenbergChainEvaluator(HeisenbergChainParameters& p){
  
  this->p = p;
}

HeisenbergChainEvaluator::~HeisenbergChainEvaluator(){

}

Eigen::MatrixXd HeisenbergChainEvaluator::get_hamiltonian_dense(vector<SpinState>& basis){
  
  const int nstates = basis.size();
  Eigen::MatrixXd h(nstates, nstates);
  h = Eigen::MatrixXd::Zero(nstates, nstates);
  fptype element=0;
  #pragma omp parallel for schedule(dynamic,1)
  for(int i=0;i<nstates;i++){
    for(int j=i;j<nstates;j++){ //only calculate upper half of the matrix, because it is hermitian
      element = this->get_element(&basis[i], &basis[j]);
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

Eigen::SparseMatrix<fptype> HeisenbergChainEvaluator::get_hamiltonian_sparse(vector<SpinState>& basis, fptype element_threshold){
  
  const int nstates = basis.size();
  Eigen::SparseMatrix<fptype> h(nstates, nstates);
  fptype element=0;
  for(int i=0;i<nstates;i++){
    for(int j=i;j<nstates;j++){ //only calculate upper half of the matrix, because it is hermitian
      element = this->get_element(&basis[i], &basis[j]);
      if(fabs(element) > element_threshold){ //make matrix sparse by checking whether element magnitude exceeds acertain threshold
        if(i == j){
	  h.insert(i,i) = element;
        }
        else{
	  h.insert(i,j) = element;
	  h.insert(j,i) = element;
	}
      }
    }
  }
  
  h.makeCompressed();
  
  return h;
}

fptype HeisenbergChainEvaluator::get_element(SpinState* u, SpinState* v){
  
  const int Nsites = u->statevec.size();
  fptype element = 0;
  for(int k=1;k<Nsites;k++){
    element += -p.J/2.0*(v->dot(u->apply_sminus(k)->apply_splus(k-1)) + v->dot(u->apply_sminus(k-1)->apply_splus(k))) - p.J*v->dot(u->apply_sz(k-1)->apply_sz(k));
  }
  if(p.periodic){
    element += -p.J/2.0*(v->dot(u->apply_sminus(0)->apply_splus(Nsites-1)) + v->dot(u->apply_sminus(Nsites-1)->apply_splus(0))) - p.J*v->dot(u->apply_sz(Nsites-1)->apply_sz(0));
  }
  
  //add magnetic field term
  for(int k=0;k<Nsites;k++){
    element -= p.g*p.B*v->dot(u->apply_sz(k));
  }
  
  return element;
};

LanczosSolver::LanczosSolver(){
  
}

LanczosSolver::~LanczosSolver(){
  
}

void LanczosSolver::compute(Eigen::SparseMatrix<fptype> m, int nev, const fptype ev_threshold){
  
  const int maxsteps = 50;
  const int matdim = m.cols();
  nev = min(nev, matdim-1);
  hessenberg = Eigen::MatrixXd::Zero(maxsteps, maxsteps);
  eigenvalues = Eigen::VectorXd::Zero(nev);
  
  Eigen::VectorXd lv_laststep, lv_lastlaststep;
  Eigen::VectorXd lv = Eigen::VectorXd::Random(matdim); //start iteration vector randomly and normalize it
  lv.normalize();
  
  hessenberg(0,0) = lv.transpose() * m * lv; //evaluate d0
  lv_laststep = lv;
  lv = m * lv - hessenberg(0,0) * lv; //obtain Phi1 without normalization
  hessenberg(1,0) = lv.norm(); //evaluate f1
  lv.normalize(); 
  
  //start Lanczos iteration
  int step = 1;
  Eigen::VectorXd eigenvalues_laststep = Eigen::VectorXd::Zero(nev);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
  do{
    hessenberg(step, step) = lv.transpose() * m * lv;
    hessenberg(step-1,step) = lv_laststep.transpose() * m * lv;
    lv_lastlaststep = lv_laststep;
    lv_laststep = lv;
    lv = m * lv_laststep - hessenberg(step,step) * lv_laststep - hessenberg(step-1, step) * lv_lastlaststep;
    hessenberg(step+1,step) = lv.norm(); //this norm is a criterion for the Lanczos algorithm to stop, if it is below threshold Krylov subspace is exhausted, but when checking mind that there is a step++ statement after this
    //cout << "Norm: " << hessenberg(step+1,step) << endl;
    lv.normalize();
    
    /*if(step<5)
      cout << hessenberg.topLeftCorner(step, step) << endl << endl;*/
    
    step++;
    solver.compute(hessenberg.topLeftCorner(step, step));
    if(step > nev){
      eigenvalues_laststep = eigenvalues.head(nev);
      eigenvalues = solver.eigenvalues().head(nev);
    }
  }
  while((((eigenvalues - eigenvalues_laststep).norm() > ev_threshold) && (hessenberg(step,step-1) > 1e-7) && (step < (maxsteps-1))) || (step < (nev+2)));
  
  if((maxsteps-1) == step){
    cout << "Lanczos algorithm did not converge." << endl;
  }
  else{
    cout << "Lanczos algorithm converged after " << step << " steps." << endl;
  }
}

vector<int> get_allowed_sz(const int nsites){
  //sz is the net spin measured in units of 1/2
  vector<int> sz;
  for(int i=-nsites;i<(nsites+1);i+=2){
    sz.push_back(i);
  }
  
  return sz;
}