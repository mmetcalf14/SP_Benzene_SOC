//
//  Hamiltonian.cpp
//  ActandAro_SOC_SPM
//
//  Created by mekena McGrew on 1/9/17.
//  Copyright © 2017 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(double _t1, double _t2, double _gamma_SOC, int Nsites)
{
    t1 = _t1;
    t2 = _t2;
    gamma_SOC = _gamma_SOC;
    L = Nsites;
    
    Ham_Mat = MatrixXcd::Zero(2*L, 2*L);
    Ham_Top = MatrixXcd::Zero(L,L);
    Ham_Btm = MatrixXcd::Zero(L,L);
}

MatrixXcd Hamiltonian::Block_DenseHam_Periodic(int species, double phi)
{
    MatrixXcd H = MatrixXcd::Zero(L,L);
    
    for(int i = 0; i < L; i++)
    {
        
        if( (i%2) == 0)
        {
            
            if( (i+1) < L)
            {
                H(i+1,i) = -t1;
                H(i,i+1) = -t1;
            }
            else if ((i+1 == L))
            {

                H(0,i) = -t1;
                H(i,0) = -t1;
            }
            else{
                
            }
            
        }
        else
        {
            if( (i+1) < L)
            {
                
                H(i+1,i) = -t2;
                H(i,i+1) = -t2;
            }
            else if ((i+1 == L))
            {
                H(0,i) = -t2;
                H(i,0) = -t2;
            }
            else{
                
            }
            
        }
    }
    
    complex<double> val;

    //spin-dependent SOC term
    if(species == 0)//spin up, top block
    {
        val = gamma_SOC*exp(I*phi);
//            if(val.real() < 1e-13)
//            {
//                val.real(0);
//            }
        //H(0,2) = I*gamma_SOC;
        //H(2,0) = -1.*I*gamma_SOC;
        H(2,0) = val;
        
        H(0,2) = conj(val);
    }
    else if(species == 1)//spin down, bottom block
    {
        val = -1.*gamma_SOC*exp(I*phi);
//            if(val.real() < 1e-13)
//            {
//                val.real(0);
//            }
//        H(0,2) = -1.*I*gamma_SOC;
//        H(2,0) =  I*gamma_SOC;
        H(2,0) = val;
        H(0,2) = conj(val);
    }
    else{
        cout << "Wrong spin species \n";
    }

    //cout << "Block Ham:\n" << H << endl;
    return H;
}

MatrixXcd Hamiltonian::BlockDense_TotalSOC(int species, double phi)
{
    MatrixXcd H = MatrixXcd::Zero(L,L);
    complex<double> val;
    if(species == 0)
    {
        val = gamma_SOC*exp(I*phi);
    }
    else{
        val = -1.*gamma_SOC*exp(I*phi);
    }
    
    for(int i = 0; i < L; i++)
    {
            if( (i+1) < L)//taking care of NN hopping
            {
                H(i+1,i) = -t1;
                H(i,i+1) = -t1;
            
            }
            else if ((i+1 == L))
            {
                
                H(0,i) = -t1;
                H(i,0) = -t1;
            }
            else{
                
            }
        
        if( (i+2) < L)//taking care of NNN hopping
        {
            H(i+2,i) = conj(val);
            H(i,i+2) = val;
            //cout << "i: " << i << " i+2: " << i+2 << endl;
            
        }
        else if ((i+2 == L))
        {
            //cout << "ii: " << i << " io+2: " << i+2 << endl;
            H(0,i) = conj(val);
            H(i,0) = val;
        }
        else{
            //cout << "i0: " << i << " ii+2: " << i+2 << endl;
            H(1,i) = conj(val);
            H(i,1) = val;
        }
    }

    return H;
  
}
void Hamiltonian::Peierls_Hamiltonian_pb()
{
    std::vector<Tp> TL;
    complex<double> val = -t1*exp(I*Pi/4.);
    
    cout << "Building Triplets\n";
    for(int i = 0; i < L; i++)
    {
        cout << "i: " << i << endl;
        if( (i+1) < L)
        {
            TL.push_back(Tp(i+1,i,val));//J1 == 1
            TL.push_back(Tp(i,i+1,conj(val)));
        }
        else if ((i+1 == L))
        {
            //                Ham_Mat(0, i) = -J1;
            //                Ham_Mat(i, 0) = -J1;
            TL.push_back(Tp(0,i,val));
            TL.push_back(Tp(i,0,conj(val)));
        }
        else{
            
        }
        
        
    }
    cout << "Setting Ham\n";
    Ham_Sp.setFromTriplets(TL.begin(), TL.end());
    //cout << "Hamiltonian: \n" << Ham_Mat_C << endl;
    cout << "Diagonalizing\n";
    Diagonalize_CompHam_SP(Ham_Sp);
    
    
}


SpMat Hamiltonian::Block_SpHam_Periodic(int species)
{
    //one option is to build the top block and the bottom block then construct Ham
    std::vector<Tp> TL;
    SpMat H;
    H.resize(L,L);
    
    //Hopping along ring
    for(int i = 0; i < L; i++)
    {
        
        if( (i%2) == 0)
        {
            
            if( (i+1) < L)
            {
                TL.push_back(Tp(i+1,i,-t1));
                TL.push_back(Tp(i,i+1,-t1));
            }
            else if ((i+1 == L))
            {
                //                Ham_Mat(0, i) = -J1;
                //                Ham_Mat(i, 0) = -J1;
                TL.push_back(Tp(0,i,-t1));
                TL.push_back(Tp(i,0,-t1));
            }
            else{
                
            }
            
        }
        else
        {
            if( (i+1) < L)
            {
                
                TL.push_back(Tp(i+1,i,-t2));
                TL.push_back(Tp(i,i+1,-t2));
            }
            else if ((i+1 == L))
            {
                TL.push_back(Tp(0,i,-t2));
                TL.push_back(Tp(i,0,-t2));
            }
            else{
                
            }
            
        }
    }
    
    
    //spin-dependent SOC term
    //This Hamiltonian is missing an onsite term and a Rashba term if compared directly to the rice mele
    //in the Rice Mele model gamma = 0.06t
    if(species == 0)//spin up, top block
    {
        TL.push_back(Tp(0,4,I*gamma_SOC));
        TL.push_back(Tp(4,0,-1.*I*gamma_SOC));
    }
    else if(species == 1)//spin down, bottom block
    {
        TL.push_back(Tp(0,4,-1.*I*gamma_SOC));
        TL.push_back(Tp(4,0, I*gamma_SOC));
    }
    else{
        cout << "Wrong spin species \n";
    }
    
    
    H.setFromTriplets(TL.begin(), TL.end());
    
    //cout << "Block Ham: \n" << H << endl;

    return H;
}

void Hamiltonian::Total_Ham(int Ltp, double phi)
{
 
    if(Ltp == 0)
    {
    Ham_Top = Block_DenseHam_Periodic(0, phi);//up
    Ham_Btm = Block_DenseHam_Periodic(1, phi);//dn
    }
    else if(Ltp == 1)//in the future LTp could represent the # of links
    {
        Ham_Top = BlockDense_TotalSOC(0, phi);
        Ham_Btm = BlockDense_TotalSOC(1, phi);
    }
    else{
        
    }
    
    Ham_Mat.topLeftCorner(L,L) = Ham_Top;
    Ham_Mat.bottomRightCorner(L,L) = Ham_Btm;
//    Ham_Top.setZero();
//    Ham_Btm.setZero();
    
    //cout << "Total Ham \n" <<Ham_Mat << endl;
//    cout << "Ham top\n" <<Ham_Top << endl;
//    cout << "Ham btm\n" <<Ham_Btm << endl;
    
//    Ham_Sp = Ham_Mat.sparseView();
//    
//    cout << "Sparse Ham \n" <<Ham_Sp << endl;
    //Ham_Mat.setZero();
    
    //Do not do block diagonalization. Want all eigenvectors. Although you could...
    //Diagonalize_CompHam(Ham_Mat);
    Diagonalize_CompHam(Ham_Top);
}

void Hamiltonian::Diagonalize(ofstream &output1, ofstream &output2)
{
    MatrixXcd EVec_up;
    VectorXd EVal_up;
    MatrixXcd EVec_dn;
    VectorXd EVal_dn;
    
    //SelfAdjointEigenSolver<Eigen::MatrixXcd> Diag(Ham_Mat);
    SelfAdjointEigenSolver<Eigen::MatrixXcd> Diag_up(Ham_Top);//spin up
    SelfAdjointEigenSolver<Eigen::MatrixXcd> Diag_dn(Ham_Btm);//spin down
    //ComplexEigenSolver<Eigen::MatrixXcd> Diag(Ham_Mat);
    EVal_up = Diag_up.eigenvalues();//eigen values are real, eigen vec not
    EVec_up = Diag_up.eigenvectors();
    EVal_dn = Diag_dn.eigenvalues();//eigen values are real, eigen vec not
    
    EVec_dn = Diag_dn.eigenvectors();
    
   // cout << "Eigenvalues: " << EVal << endl;
    
//    for(int i = 0; i< 2*L; i++)
//    {
//        //output << i << "\t" << EVal(i).real() << endl;
//        output << i << "\t" << EVal(i) << endl;
//    }
    
//        for(int i = 0; i< L; i++)
//        {
//            //output << i << "\t" << EVal(i).real() << endl;
//            output << i << " " << EVal_up(i) << " " << EVal_dn(i) << endl;
//            cout << i << " " << EVal_up(i) << " " << EVal_dn(i) << endl;
//            //THese Eigenvalues pair and are the same.
//        }

    for(int i = 0; i < L; i++)
    {
        if(i == 0)
        {
            output1 << gamma_SOC << " " << EVal_up(i) << " " ;
            output2 << gamma_SOC << " " << EVal_dn(i) << " " ;
        }
        else{
            output1 << EVal_up(i) << " ";
            output2 << EVal_dn(i) << " ";
        }
    }
    
    output1 << endl;
    output2 << endl;
    
//
//    VectorXd Vec = EVal.real();
//    cout <<"Vec\n"<< Vec << endl;
}

void Hamiltonian::Diagonalize_CompHam_SP(SpMat Ham)
{
    cout << "setting ham\n";
    MatrixXcd H = MatrixXcd(Ham);
    cout << "Diagonalizing\n";
    SelfAdjointEigenSolver<Eigen::MatrixXcd> Diag(H);
    EVal = Diag.eigenvalues();
    EVec = Diag.eigenvectors();
}

void Hamiltonian::Diagonalize_CompHam(MatrixXcd H)
{
   cout << "Diagonalizing\n";
    SelfAdjointEigenSolver<Eigen::MatrixXcd> Diag(H);
    EVal = Diag.eigenvalues();
    EVec = Diag.eigenvectors();
//    cout << "Evals:\n" << EVal << endl;
//    cout << "EVec: " << EVec << endl;
    //GS = EVec.col(0);
    
    //cout << "GS: " << GS << endl;
//    cout << "GS mod: \n";
//    for(int it = 0; it < L; it++)
//    {
//        cout << GS.col(0).row(it).adjoint()*GS.col(0).row(it) << endl;
//    }
}

complex<double> Hamiltonian::BondCurrent(int type, int site_i, int site_j)
{
    complex<double> J;
    complex<double> Cij;
    
    Cij = conj(EVec(site_i, type))*EVec(site_j, type)*Ham_Mat(site_i,site_j);
    J = -2.*Cij.imag();
    
    return J;
}

VectorXd Hamiltonian::ReturnEval()
{
        return EVal;
}

void Hamiltonian::ResetGamma(double g)
{
    gamma_SOC = g;
}

void Hamiltonian::ResetHam()
{
    Ham_Mat.setZero();
    Ham_Top.setZero();
    Ham_Btm.setZero();
}





