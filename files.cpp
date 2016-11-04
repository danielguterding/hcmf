//files.cpp
#include <iostream>

#include "files.hpp"

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

HeisenbergBondReader::HeisenbergBondReader(){
  
}

void HeisenbergBondReader::read_file(const string infilename){
  
  this->bonds.resize(0);
  
  boost::filesystem::path infilepath(infilename);
  if(!boost::filesystem::exists(infilepath)){cout << "Warning! Input file does not exist!" << endl;};
  boost::filesystem::ifstream infilehandle(infilepath);
  
  string line = "";
  vector<string> splitline;
  
  getline(infilehandle, line); //discard first two lines
  getline(infilehandle, line);
  while(getline(infilehandle, line)){
    line = trim_all(line); //erase trailing and leading spaces, reduce intermediate spaces to one space
    boost::split(splitline, line, boost::is_any_of("\t "));
    if(splitline.size()<3){ //break loop after all Heisenberg bonds have been read
      break;
    }
    HeisenbergBond b;
    b.s1 = boost::lexical_cast<int>(splitline[0]);
    b.s2 = boost::lexical_cast<int>(splitline[1]);
    b.J = boost::lexical_cast<fptype>(splitline[2]);
    this->bonds.push_back(b);
  }
  infilehandle.close();
  
  sort(this->bonds.begin(), this->bonds.end(), [](const HeisenbergBond& b1, const HeisenbergBond& b2){return b1.s2 > b2.s2;}); //sort bonds after s2 with inline lambda function
  this->nsites = this->bonds[0].s2+1;
}

SpinBasisWriter::SpinBasisWriter(){
  
}

void SpinBasisWriter::write_basis(string outfilename, vector<SpinState> basis){
  
  boost::filesystem::path outfilepath(outfilename);
  boost::filesystem::ofstream outfilehandle(outfilepath);
  
  outfilehandle << "#state idx, bitpattern" << endl;
  for(uint i=0;i<basis.size();i++){
    SpinState s = basis[i];
    string outstr = "";
    for(uint j=0;j<s.statevec.size();j++){
      if(s.statevec[j]){
        outstr += "1";
      }
      else{
        outstr += "0";
      }
    }
    outfilehandle << boost::format("% 14i ") % i << outstr << endl;
  }
  
  outfilehandle.close();
}

string trim_all(const std::string &str){  //with a more recent version of boost boost::trim_all() can be used instead of this function

  return boost::algorithm::find_format_all_copy(
    boost::trim_copy(str),
    boost::algorithm::token_finder (boost::is_space(),boost::algorithm::token_compress_on),
    boost::algorithm::const_formatter(" "));
}