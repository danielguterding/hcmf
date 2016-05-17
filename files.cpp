//files.cpp

#include "files.hpp"

HeisenbergBondReader::HeisenbergBondReader(){
  
}

void HeisenbergBondReader::read_file(string infilename){
  
  this->bonds.resize(0);
  
  boost::filesystem::path infilepath(infilename);
  if(!boost::filesystem::exists(infilepath)){cout << "Warning! Input file does not exist!" << endl;};
  boost::filesystem::ifstream infilehandle(infilepath);
  
  string line = "";
  vector<string> splitline;
  
  getline(infilehandle, line); //discard first line
  while(getline(infilehandle, line)){
    line = trim_all(line); //erase trailing and leading spaces, reduce intermediate spaces to one space
    boost::split(splitline, line, boost::is_any_of("\t "));
    HeisenbergBond b;
    b.s1 = boost::lexical_cast<int>(splitline[0]);
    b.s2 = boost::lexical_cast<int>(splitline[1]);
    b.J = boost::lexical_cast<fptype>(splitline[2]);
    this->bonds.push_back(b);
  }
  infilehandle.close();
}

SiteDependentMagneticFieldReader::SiteDependentMagneticFieldReader(){
  
}

void SiteDependentMagneticFieldReader::read_file(string infilename){
  
  this->fields.resize(0);
  
  boost::filesystem::path infilepath(infilename);
  if(!boost::filesystem::exists(infilepath)){cout << "Warning! Input file does not exist!" << endl;};
  boost::filesystem::ifstream infilehandle(infilepath);
  
  string line = "";
  vector<string> splitline;
  
  getline(infilehandle, line); //discard first line
  while(getline(infilehandle, line)){
    line = trim_all(line); //erase trailing and leading spaces, reduce intermediate spaces to one space
    boost::split(splitline, line, boost::is_any_of("\t "));
    SiteDependentMagneticField f;
    f.s = boost::lexical_cast<int>(splitline[0]);
    f.h = boost::lexical_cast<fptype>(splitline[1]);
    this->fields.push_back(f);
  }
  infilehandle.close();
}

string trim_all(const std::string &str){  //with a more recent version of boost boost::trim_all() can be used instead of this function

  return boost::algorithm::find_format_all_copy(
    boost::trim_copy(str),
    boost::algorithm::token_finder (boost::is_space(),boost::algorithm::token_compress_on),
    boost::algorithm::const_formatter(" "));
}