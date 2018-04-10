 
/*
 * hs071_nlp.hpp
 *
 *  Created on: Dec 3, 2017
 *      Author: martin
 */

#ifndef PROTDATA_HPP_
#define PROTDATA_HPP_


#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <ctime>


using namespace std;

class ProtData {
public:
    
    ProtData();
    ProtData(const string charge, const string dielx, const string diely, const string dielz, const string kappa, const string pot);

    void extract_delta(string& str, double &d1, double &d2, double &d3);
    void extract_orig(string& str, double &o1, double &o2, double &o3);
    void extract_dim(string& str, unsigned &n1, unsigned &n2, unsigned &n3);
    void read_dxfile(const string &file, unsigned& n1, unsigned &n2, unsigned &n3, double &o1, double &o2, double &o3, double &d1, double &d2, double &d3, vector<vector<vector<double>>>  &v);
    void freeMemory();
    void write_file(const string filename, vector<vector<vector<double>>> &v);
    void fix_vector(vector<vector<vector<double>>> &vec, const unsigned type);


    /** Default destructor */
    
	//~ProtData();

	/*
	 * nx, ny, nz dimensions
	 * ox. oy. iz coordinates of origin
	 * hx, hy, hz  grid widths
	 */
	
	unsigned nx, ny, nz, N;
	double ox, oy, oz;
	double hx, hy, hz, h;
        vector<vector<vector<double>>> charge, dielx, diely, dielz, kappa, pot;
        
        vector<unsigned> rows, cols;
        vector<float> vals;
        
        
	
private:


};



#endif /* PROTDATA_HPP_ */
