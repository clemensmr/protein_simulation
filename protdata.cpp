 
/*
 * matrixop.cpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#include "protdata.hpp"

using namespace std;

// constructor

PROTDATA::PROTDATA() {
}

PROTDATA::PROTDATA(string charge_file, string dielx_file, string diely_file, string dielz_file, string kappa_file, string pot_file){
    


}

const string commentary("#");
const string delim(" ,");

void PROTDATA::extract_delta(string& str, double &d1, double &d2, double &d3)
{
  // get rid of delimiters 
  size_t found = str.find_first_not_of(delim);
  if(found != string::npos) str.erase(0, found);

  // extract first dimension
  found = str.find_first_of(delim);
  const string value1 = str.substr(0, found);
  str.erase(0, value1.size());
  double t = stof(value1);
  if (t!=0.0) d1 = t;

   // get rid of delimiters 
  found = str.find_first_not_of(delim);
  if(found != string::npos) str.erase(0, found);

  // extract second dimension
  found = str.find_first_of(delim);
  const string value2 = str.substr(0, found);
  str.erase(0, value2.size());
  t = stof(value2);
  if (t!=0.0) d2 = t;

  // get rid of delimiters 
  found = str.find_first_not_of(delim);
  if(found != string::npos) str.erase(0, found);

  // extract third dimension
  found = str.find_first_of(delim);
  const string value3 = str.substr(0, found);
  str.erase(0, value3.size());
  t = stof(value3);
  if (t!=0.0) d3 = t;
}



void PROTDATA::extract_orig(string& str, double &o1, double &o2, double &o3)
{
  // get rid of delimiters 
  size_t found = str.find_first_not_of(delim);
  if(found != string::npos) str.erase(0, found);

  // extract first dimension
  found = str.find_first_of(delim);
  const string value1 = str.substr(0, found);
  str.erase(0, value1.size());
  o1 = stof(value1);

   // get rid of delimiters 
  found = str.find_first_not_of(delim);
  if(found != string::npos) str.erase(0, found);

  // extract second dimension
  found = str.find_first_of(delim);
  const string value2 = str.substr(0, found);
  str.erase(0, value2.size());
  o2 = stof(value2);

  // get rid of delimiters 
  found = str.find_first_not_of(delim);
  if(found != string::npos) str.erase(0, found);

  // extract third dimension
  found = str.find_first_of(delim);
  const string value3 = str.substr(0, found);
  str.erase(0, value3.size());
  o3 = stof(value3);

  //cout << "o1=" << o1 << ", o2="<< o2 << ", o3=" << o3 << endl;
}


void PROTDATA::extract_dim(string& str, unsigned &n1, unsigned &n2, unsigned &n3)
{
  // get rid of delimiters 
  size_t found = str.find_first_not_of(delim);
  if(found != string::npos) str.erase(0, found);

  // extract first dimension
  found = str.find_first_of(delim);
  const string value1 = str.substr(0, found);
  str.erase(0, value1.size());
  n1 = stoi(value1);

   // get rid of delimiters 
  found = str.find_first_not_of(delim);
  if(found != string::npos) str.erase(0, found);

  // extract second dimension
  found = str.find_first_of(delim);
  const string value2 = str.substr(0, found);
  str.erase(0, value2.size());
  n2 = stoi(value2);

  // get rid of delimiters 
  found = str.find_first_not_of(delim);
  if(found != string::npos) str.erase(0, found);

  // extract third dimension
  found = str.find_first_of(delim);
  const string value3 = str.substr(0, found);
  str.erase(0, value3.size());
  n3 = stoi(value3);

  cout << "n1=" << n1 << ", n2="<< n2 << ", n3=" << n3 << endl;
}


void PROTDATA::read_dxfile(const string &file, unsigned& n1, unsigned &n2, unsigned &n3, double &o1, double &o2, double &o3, double &d1, double &d2, double &d3, vector<double>  &v)
{
  // open file and check if it was successful

  ifstream infile(file.data(), ifstream::in);
  if(!infile.is_open()) {
    cerr << file.data() << ": error: could not open file!" << endl;  
    exit(1);
  }

  cout << "Reading " << file.data() << endl;
  
  string line; // containing read line

  // read file line by line until end of file is reached
  while(!infile.eof()) {
        getline(infile, line);

        // get rid of comments
        size_t found = line.find_first_of(commentary);
        if(found != string::npos) line.erase(found, string::npos);
        //    cout << "com: " << line << endl;

        // extract dimensions
        const string s1("object 1 class gridpositions counts");
        found = line.find(s1);
        if (found != string::npos) {
        line.erase(0, s1.size());
        extract_dim(line, n1, n2, n3);
        }
        const string s2("origin");
        found = line.find(s2);
        if (found != string::npos) {
        line.erase(0, s2.size());
        extract_orig(line, o1, o2, o3);
        cout << "o1=" << o1 << ", o2="<< o2 << ", o3=" << o3 << endl;
        }
        const string s3("delta");
        found = line.find(s3);
        if (found != string::npos) {
        line.erase(0, s3.size());
        extract_delta(line, d1, d2, d3);
        if (d1!=0.0 && d2!=0.0 && d3!=0.0)  cout << "d1=" << d1 << ", d2="<< d2 << ", d3=" << d3 << endl;
        }

        const string s4("object 3 class array type");
        found = line.find(s4);
        if (found != string::npos) {
        //      cout << line << endl;
        const unsigned N = n1*n2*n3;

        v.resize(N);
        unsigned i;
        for (i=0; i<N; i+=3) {
        getline(infile, line);
        //        cout << line << endl;
        extract_orig(line, v[i], v[i+1], v[i+2]);
        }
        cout << "Read " << i << " values of " << N << "." << endl;
        
    }
    
  }

  // close file
  infile.close();
    
}