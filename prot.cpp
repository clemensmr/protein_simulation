 
 /*
  * main.cpp
  *
  *  Created on: Dec 3, 2017
  *      Author: martin
  */
 #include "protdata.hpp"
 #include "CSmatrix.h"
 #include "args.hxx"
 #include "CS_helper.h"
 #include "solvers.h"
 
 #include <string>
 #include <vector>
 #include <time.h>
 #include <fstream>
 #include <cmath>
 
 
 
 using namespace AHMED;
 
 void fillMatrix(std::vector<AccumulatingMap<float>> &csMatrix, const ProtData &data);
 void fill_rhs(float *&rhs, const ProtData &data);
 void writeMatrix(const std::vector<AccumulatingMap<double>> csMatrix);
 void writeCRSMatrix(std::string filename, const unsigned n,const unsigned* iA, const unsigned* jA, const float* A);
 void writeArray(const double *vec, const unsigned n);
 
 
 int main(int argc, char** argv) {
     int MPC_horizon = 100;
     std::string CRSname = "CRS_test";
     
     string pot = "/home/clemens/data/APBS/cytc-pot-0.12-1.0.dx";
     string kappa = "/home/clemens/data/APBS/cytc-kappa.dx.dx";
     string dielx = "/home/clemens/data/APBS/cytc-dielx.dx.dx";
     string diely = "/home/clemens/data/APBS/cytc-diely.dx.dx";
     string dielz = "/home/clemens/data/APBS/cytc-dielz.dx.dx";
     string charge = "/home/clemens/data/APBS/cytc-charge.dx.dx";
     
     args::ArgumentParser parser("convection diffusion equation 1d.", "This goes after the options.");
     args::HelpFlag help(parser, "help", "Display this help menu",{'h', "help"});
     args::CompletionFlag completion(parser,{"complete"});
     
     args::Flag test(parser, "test", "toggle test files",{'t', "test"});
     
     args::ValueFlag<int> MPC_horizon_(parser, "mpc horizon", "MPC horizon, N >= 1",{'N', "mpc"});
     
     
     try {
         parser.ParseCLI(argc, argv);
     }
     catch (args::Completion e) {
         std::cout << e.what();
         return 0;
     }
     catch (args::Help) {
         std::cout << parser;
         return 0;
     }
     catch (args::ParseError e) {
         std::cerr << e.what() << std::endl;
         std::cerr << parser;
         return 1;
     }
     
     //int
     if (test) {
         pot = "/home/clemens/data/APBS/pot_test";
         kappa = "/home/clemens/data/APBS/kappa_test";
         dielx = "/home/clemens/data/APBS/dielx_test";
         diely = "/home/clemens/data/APBS/diely_test";
         dielz = "/home/clemens/data/APBS/dielz_test";
         charge = "/home/clemens/data/APBS/charge_test";
     }
     
     std::cout << "this works" << std::endl;
     
     unsigned N;
     float *A, *rhs;
     unsigned *iA, *jA;
     
     try{
         ProtData data(charge, dielx, diely, dielz, kappa, pot);
         
         N = data.N;
         
         std::vector<AccumulatingMap<float>> csMatrix(data.N);
         
         fillMatrix(csMatrix, data);
         std::cout << "matrix filled" << std::endl;
         fill_rhs(rhs, data);
         std::cout << "right side filled" << std::endl;
         
         data.freeMemory();
         //std::cout << "memory free" << std::endl;
         
         convertToCRSMatrix(csMatrix, iA, jA, A);
         std::cout << "crs matrix created" << std::endl;
         
         
     }
     catch(const std::exception& e){
     }
     
     //writeCRSMatrix(CRSname, N, iA, jA, A);
     
     AHMED::CRSMatrix<float> mat(N);
     
     mat.noRows = N;
     mat.iA = iA;
     mat.jA = jA;
     mat.A = A;
     //solve system
     std::pair<float, unsigned> re;
     
     //re = CG(const Matrix<T>& A, const typename num_traits<T>::abs_type eps, const unsigned nsteps, const T * const b, T * const x);
     //writeMatrix(csMatrix);
     
     
     //writeArray(rhs, data.N);
     
     delete[] A;
     delete[] iA;
     delete[] jA;
     delete[] rhs;
     
     return 0;
 }
 
 
 void fillMatrix(std::vector<AccumulatingMap<float>> &csMatrix,  const ProtData &data){
     
     
     double nx = data.nx;
     double ny = data.ny;
     double nz = data.nz;
     double h = data.h;
     
     
     //diagonal    
     for(unsigned i = 0; i < nx; ++i){
         for(unsigned k = 0; k < ny; ++k){
             for(unsigned l = 0; l < nz; ++l){
                 
                 double kappa = data.kappa[i][k][l];
                 csMatrix[i*ny*nz+k*nz+l].accumulatedPush(kappa * kappa * h * h, i*ny*nz+k*nz+l);   
                 
                 if(i != 0){
                     csMatrix[i*ny*nz+k*nz+l].accumulatedPush(data.dielx[i][k][l], i*ny*nz+k*nz+l);
                 }
                 if(i != nx-1){
                     csMatrix[i*ny*nz+k*nz+l].accumulatedPush(data.dielx[i+1][k][l], i*ny*nz+k*nz+l);
                 }
                 if(k != 0){
                     csMatrix[i*ny*nz+k*nz+l].accumulatedPush(data.diely[i][k][l], i*ny*nz+k*nz+l);
                 }
                 if(k != ny-1){
                     csMatrix[i*ny*nz+k*nz+l].accumulatedPush(data.diely[i][k+1][l], i*ny*nz+k*nz+l);
                 }
                 if(l != 0){
                     csMatrix[i*ny*nz+k*nz+l].accumulatedPush(data.dielz[i][k][l], i*ny*nz+k*nz+l);
                 }
                 if(l != nz-1){
                     csMatrix[i*ny*nz+k*nz+l].accumulatedPush(data.dielz[i][k][l+1], i*ny*nz+k*nz+l);
                 }
                 
             }
         }
     }
     
     //outer band -> dielx
     for(unsigned i = 1; i < nx; ++i){
         for(unsigned k = 0; k < ny; ++k){
             for(unsigned l = 0; l < nz; ++l){
                 //lower
                 csMatrix[i*ny*nz+k*nz+l].accumulatedPush(-data.dielx[i][k][l], (i-1)*ny*nz+k*nz+l);
                 //upper
                 csMatrix[(i-1)*ny*nz+k*nz+l].accumulatedPush(-data.dielx[i][k][l], i*ny*nz+k*nz+l);
             }
         }
     }
     
     //middle band -> diely
     for(unsigned i = 0; i < nx; ++i){
         for(unsigned k = 1; k < ny; ++k){
             for(unsigned l = 0; l < nz; ++l){
                 //lower
                 csMatrix[i*ny*nz+k*nz+l].accumulatedPush(-data.diely[i][k][l], i*ny*nz+(k-1)*nz+l);
                 //upper
                 csMatrix[i*ny*nz+(k-1)*nz+l].accumulatedPush(-data.diely[i][k][l], i*ny*nz+k*nz+l);
             }
         }
     }
     
     //inner band -> dielz
     for(unsigned i = 0; i < nx; ++i){
         for(unsigned k = 0; k < ny; ++k){
             for(unsigned l = 1; l < nz; ++l){
                 //lower
                 csMatrix[i*ny*nz+k*nz+l].accumulatedPush(-data.dielz[i][k][l], i*ny*nz+k*nz+(l-1));
                 //upper
                 csMatrix[i*ny*nz+k*nz+(l-1)].accumulatedPush(-data.dielz[i][k][l], i*ny*nz+k*nz+l);
             }
         }
     }
 }
 
 void fill_rhs(float *&rhs, const ProtData &data){
     rhs = new float[data.N];
     
     unsigned nx = data.nx;
     unsigned ny = data.ny;
     unsigned nz = data.nz;
     
     for(unsigned i = 0; i < nx; ++i){
         for(unsigned k = 0; k < ny; ++k){
             for(unsigned l = 0; l < nz; ++l){
                 rhs[i*ny*nz+k*ny+l] = 4.*M_PI/data.h*data.charge[i][k][l];
                 
                 if(i == 0){
                     rhs[i*ny*nz+k*ny+l] += data.dielx[i][k][l] * data.pot[i][k+1][l+1];
                 }
                 if(i == nx-1){
                     rhs[i*ny*nz+k*ny+l] += data.dielx[i+1][k][l] * data.pot[i+2][k+1][l+1];
                 }
                 if(k == 0){
                     rhs[i*ny*nz+k*ny+l] += data.diely[i][k][l] * data.pot[i+1][k][l+1];
                 }
                 if(k == ny-1){
                     rhs[i*ny*nz+k*ny+l] += data.diely[i][k+1][l] * data.pot[i+1][k+2][l+1];
                 }
                 if(l == 0){
                     rhs[i*ny*nz+k*ny+l] += data.dielz[i][k][l] * data.pot[i+1][k+1][l];
                 }
                 if(l == nz-1){
                     rhs[i*ny*nz+k*ny+l] += data.dielz[i][k][l+1] * data.pot[i+1][k+1][l+2];
                 }
             }
         }
     }
 }
 
 
 void writeMatrix(const std::vector< AccumulatingMap<double>> csMatrix){
     std::ofstream out("../../testmatrix");
     
     if(!out.is_open()){
         std::cout << "couldn't open output file" << std::endl;
         exit(1);
     }
     
     for(unsigned i = 0; i < csMatrix.size(); ++i){
         //for(unsigned k = 0; k <csMatrix[i].size(); ++k){
         for(auto it : csMatrix[i]){
             out << it.first << ", " << it.second << "\t";
         }
         
         out << std::endl;
         
     }
     out.close();
 }
 
 void writeArray(const double *vec, const unsigned n){
     std::ofstream out("../../testrhs");
     
     if(!out.is_open()){
         std::cout << "couldn't open output file" << std::endl;
         exit(1);
     }
     
     for(unsigned i = 0; i < n; ++i){
         
         out << vec[i] << std::endl;
         
     }
     out.close();
 }
 
 void writeCRSMatrix(std::string filename, const unsigned n,const unsigned* iA, const unsigned* jA, const float* A){
     std::ofstream ofs(filename);
     if(ofs.is_open()){
         for(unsigned i = 0; i <= n; ++i){
             ofs << iA[i] << " ";
         }
         for(unsigned i = 0; i < iA[n]; ++i){
             ofs << jA[i] << " ";
         }
         for(unsigned i = 0; i < iA[n]; ++i){
             ofs << A[i] << " ";
         }
         
         ofs.close();
     }
     else{
         std::cout << "couldn't open output file for CRS matrix" << std::endl;
     }
 }
 
 
 