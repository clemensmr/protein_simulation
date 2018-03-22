 
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

class PROTDATA {
public:
    
    PROTDATA();
    PROTDATA(string charge, string dielx, string diely, string dielz, string kappa, string pot);

    void build_matrix();
    void extract_delta(string& str, double &d1, double &d2, double &d3);
    void extract_orig(string& str, double &o1, double &o2, double &o3);
    void extract_dim(string& str, unsigned &n1, unsigned &n2, unsigned &n3);
    void read_dxfile(string &file, unsigned& n1, unsigned &n2, unsigned &n3, double &o1, double &o2, double &o3, double &d1, double &d2, double &d3, vector<vector<vector<double>>>  &v);
    
    void write_file(string filename, vector<vector<vector<double>>> &v);



    /** Default destructor */
    
	//~PROTDATA();

	/*
	 * nx, ny, nz dimensions
	 * ox. oy. iz coordinates of origin
	 * hx, hy, hz  grid widths
	 */
	
	unsigned nx, ny, nz;
	double ox, oy, oz;
	double hx, hy, hz;
        vector<vector<vector<double>>> charge, dielx, diely, dielz, kappa, pot;
        
        vector<unsigned> rows, cols;
        vector<double> vals;
        
        
	
private:


};



#endif /* PROTDATA_HPP_ */
