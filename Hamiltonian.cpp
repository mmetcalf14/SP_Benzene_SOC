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

MatrixXcd Hamiltonian::Block_DenseHam_Periodic(int species)
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
                //                Ham_Mat(0, i) = -J1;
                //                Ham_Mat(i, 0) = -J1;
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
    
    
    //spin-dependent SOC term
    if(species == 0)//spin up, top block
    {
        H(0,2) = I*gamma_SOC;
        H(2,0) = -1.*I*gamma_SOC;
    }
    else if(species == 1)//spin down, bottom block
    {
        H(0,2) = -1.*I*gamma_SOC;
        H(2,0) =  I*gamma_SOC;
    }
    else{
        cout << "Wrong spin species \n";
    }

    //cout << "Block Ham:\n" << H << endl;
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

void Hamiltonian::Total_Ham()
{
 
    Ham_Top = Block_DenseHam_Periodic(0);
    Ham_Btm = Block_DenseHam_Periodic(1);
    
    Ham_Mat.topLeftCorner(L,L) = Ham_Top;
    Ham_Mat.bottomRightCorner(L,L) = Ham_Btm;
//    Ham_Top.setZero();
//    Ham_Btm.setZero();
    
    //cout << "Total Ham \n" <<Ham_Mat << endl;
    
//    Ham_Sp = Ham_Mat.sparseView();
//    
//    cout << "Sparse Ham \n" <<Ham_Sp << endl;
    //Ham_Mat.setZero();
    
    //Do not do block diagonalization. Want all eigenvectors. Although you could...
    Diagonalize_CompHam(Ham_Mat);
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
    MatrixXcd H = MatrixXcd(Ham);
    SelfAdjointEigenSolver<Eigen::MatrixXcd> Diag(H);
    EVal = Diag.eigenvalues();
    EVec = Diag.eigenvectors();
}

void Hamiltonian::Diagonalize_CompHam(MatrixXcd H)
{
    SelfAdjointEigenSolver<Eigen::MatrixXcd> Diag(H);
    EVal = Diag.eigenvalues();
    EVec = Diag.eigenvectors();
    //cout << "EVec: " << EVec << endl;
    GS = EVec.col(0);
    
    cout << "GS: " << GS << endl;
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
    
    cout << "site1: " << site_i <<" site2: " << site_j << endl;
    
    //i < j
    //Cij = EVec(0, site_i)*EVec(0, site_j)*Ham_Sp.coeffRef(site_i,site_j);//sparse
    //Cij = GS(site_i)*GS(site_j)*Ham_Mat(site_i,site_j);//dense

    
    if(type == 0) //spin up
    {
        Cij = GS(site_i)*GS(site_j)*Ham_Mat(site_i,site_j);
        cout << Ham_Mat(site_i,site_j) << endl;
        cout << "Cij: " << Cij << endl;
        J = -2.*Cij.imag();
        cout << "J: " << J << endl;
    }
    else if(type == 1)//spin down
    {
        Cij = EVec(site_i, 1)*EVec(site_j, 1)*Ham_Mat(site_i,site_j);
        J = -2.*Cij.imag();
        cout << "Cij: " << Cij << endl;
        cout << "J: " << J << endl;
    }
    else{
        cout <<"somthing else\n";
    }
    
    
    return J;
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





