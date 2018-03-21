 
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

    void extract_delta(string& str, double &d1, double &d2, double &d3);
    void extract_orig(string& str, double &o1, double &o2, double &o3);
    void extract_dim(string& str, unsigned &n1, unsigned &n2, unsigned &n3);
    void read_dxfile(const string &file, unsigned& n1, unsigned &n2, unsigned &n3, double &o1, double &o2, double &o3, double &d1, double &d2, double &d3, vector<double>  &v);
    
    
    PROTDATA();
    PROTDATA(string charge, string dielx, string diely, string dielz, string kappa, string pot);



    /** Default destructor */
    ~PROTDATA();

    vector<vector<vector<double> > > charge, dielx, diely, dielz, kappa, pot;
private:


};



#endif /* PROTDATA_HPP_ */
