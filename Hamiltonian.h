//
//  Hamiltonian.h
//  ActandAro_SOC_SPM
//
//  Created by mekena McGrew on 1/9/17.
//  Copyright © 2017 Mekena Metcalf. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include "/usr/local/include/Eigen/Eigen"
#include "/usr/local/include/Eigen/Sparse"
#include "/usr/local/include/Eigen/Eigenvalues"
#include </usr/local/include/Eigen/Core>
using namespace Eigen;
using namespace std;

#ifndef Hamiltonian_h
#define Hamiltonian_h

typedef Eigen::SparseMatrix<complex<double>> SpMat;
typedef Eigen::Triplet<complex<double>> Tp;

class Hamiltonian
{
    friend class SlaterDet;
    
private:
    
    double t1, t2, gamma_SOC;
    int L;
    SpMat Ham_Sp;
//    SpMat Ham_Top;
//    SpMat Ham_Btm;
    
    MatrixXcd Ham_Mat;
    MatrixXcd Ham_Top;
    MatrixXcd Ham_Btm;
    
//    VectorXcd EVal;
    MatrixXcd EVec;
    VectorXd EVal;
    //MatrixXd EVec;
    VectorXcd GS;
    
protected:
    
public:
    
    Hamiltonian();
    Hamiltonian(double _t1, double _t2, double _gamma_SOC, int Nsites);
    
    
    SpMat Block_SpHam_Periodic(int species);
    MatrixXcd Block_DenseHam_Periodic(int species, double phi);
    MatrixXcd BlockDense_TotalSOC(int species, double phi);
    void Total_Ham(int Ltp, double phi);
    void Peierls_Hamiltonian_pb();
    void Diagonalize(ofstream &output1, ofstream &output2);
    void Diagonalize_CompHam_SP(SpMat Ham);
    void Diagonalize_CompHam(MatrixXcd H);
    void ResetGamma(double g);
    void ResetHam();
    VectorXd ReturnEval();
    complex<double> BondCurrent(int type, int site_i, int site_j);
    
};

//constants
const std::complex<double> I(0.0,1.0);
const double Pi = 4.0*atan(1);
#endif /* Hamiltonian_h */
