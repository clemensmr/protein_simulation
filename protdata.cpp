 
 /*
  * matrixop.cpp
  *
  *  Created on: Dec 3, 2017
  *      Author: martin
  */
 
 #include "protdata.hpp"
 
 using namespace std;
 
 // constructor
 
 ProtData::ProtData() {
 }
 
 ProtData::ProtData(const string charge_file, const string dielx_file, const string diely_file, const string dielz_file, const string kappa_file, const string pot_file){
     
     read_dxfile(dielx_file, nx, ny, nz, ox, oy, oz, hx, hy, hz, dielx);
     read_dxfile(diely_file, nx, ny, nz, ox, oy, oz, hx, hy, hz, diely);
     read_dxfile(dielz_file, nx, ny, nz, ox, oy, oz, hx, hy, hz, dielz);
     read_dxfile(charge_file, nx, ny, nz, ox, oy, oz, hx, hy, hz, charge);
     read_dxfile(kappa_file, nx, ny, nz, ox, oy, oz, hx, hy, hz, kappa);
     read_dxfile(pot_file, nx, ny, nz, ox, oy, oz, hx, hy, hz, pot);
     
     nx -= 2;
     ny -= 2;
     nz -= 2;
     
     N = nx*ny*nz;
     
     if(!(hx == hy && hy == hz)){
         cout << "different step widths, exiting" << endl;
         exit(1);
     }
     
     h = hx;
     
     fix_vector(kappa, 0);
     fix_vector(charge, 0);
     fix_vector(dielx, 1);
     fix_vector(diely, 2);
     fix_vector(dielz, 3);
     
     
 }
 
 const string commentary("#");
 const string delim(" ,");
 
 
 void ProtData::extract_delta(string& str, double &d1, double &d2, double &d3) {
     
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
 
 
 
 void ProtData::extract_orig(string& str, double &o1, double &o2, double &o3) {
     
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
 
 
 void ProtData::extract_dim(string& str, unsigned &n1, unsigned &n2, unsigned &n3) {
     
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
 
 
 void ProtData::read_dxfile(const string &file, unsigned& n1, unsigned &n2, unsigned &n3, double &o1, double &o2, double &o3, double &d1, double &d2, double &d3, vector<vector<vector<double>>>  &v){
     
     
     d1 = d2 = d3 = 0;
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
             
             //resize vector
             v.resize(n1);
             for(unsigned i = 0; i < n1; ++i){
                 v[i].resize(n2);
                 for(unsigned j = 0; j < n2; ++j){
                     v[i][j].resize(n3);
                 }
             }
             
             //fill vector 
             for(unsigned i = 0; i < n1; ++i){
                 for(unsigned j = 0; j < n2; ++j){
                     for(unsigned k = 0; k < n3; ++k){
                         infile >> v[i][j][k];
                     }
                 }
             }
         }
         //cout << "Read " << << " values of " << N << "." << endl;
         
     }
     
     // close file
     infile.close();
     
 }
 
 
 void ProtData::freeMemory(){
     cout << "freed memory for the first time" << endl;
     vector<vector<vector<double>>>().swap(charge);
     vector<vector<vector<double>>>().swap(pot);
     vector<vector<vector<double>>>().swap(dielx);
     vector<vector<vector<double>>>().swap(diely);
     vector<vector<vector<double>>>().swap(dielz);
     vector<vector<vector<double>>>().swap(kappa);
     
     cout << "freed memory for the first time" << endl;
 }
 
 void ProtData::write_file(const string filename, vector<vector<vector<double>>> &v){
     ofstream out(filename);
     for(unsigned i = 0; i < nx; ++i){
         for(unsigned j = 0; j < ny; ++j){
             for(unsigned k = 0; k < nz; k+=3){
                 out << v[i][j][k] << " " << v[i][j][k+1] << " " << v[i][j][k+2] << endl ;
             }
         }
     }
     out.close();
 }
 
 
 
 void ProtData::fix_vector(vector<vector<vector<double>>> &vec, const unsigned type){
     
     if(type != 1){
         vec.erase(vec.begin());
     }
     vec.erase(vec.end()-1);
     
     for(unsigned i = 0; i < vec.size(); ++i){
         
         if(type != 2){
             vec[i].erase(vec[i].begin());
         }
         vec[i].erase(vec[i].end()-1);
     }
     
     for(unsigned i = 0; i < vec.size(); ++i){
         for(unsigned k = 0; k < vec[i].size(); ++k){
             
             if(type != 3){
                 vec[i][k].erase(vec[i][k].begin());
             }
             vec[i][k].erase(vec[i][k].end()-1);
         }
     }
 }