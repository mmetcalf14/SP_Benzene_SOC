//
//  main.cpp
//  ActandAro_SOC_SPM
//
//  Created by mekena McGrew on 1/9/17.
//  Copyright © 2017 Mekena Metcalf. All rights reserved.
//

#include <iostream>
#include "Hamiltonian.h"
#include "Slater_Determinant_SPM.h"

void ProbCurrent(Hamiltonian &h);
void EnergySpec(Hamiltonian &h, ofstream &out);
void CurrWFlux(Hamiltonian &h, ofstream &out);

int main(int argc, const char * argv[]) {
    
    
    //declare constants
    double t_1 = 1.0;
    double t_2;
    double hoprat;
    double gamma;
    double dg = .01;
    int Gm = 2/dg;
    double phi_rat;
    
    int Nsite;
    int Npart;
    char output[60];
    char otherout[60];
    char currentout[80];
    //call config file eventually
    ifstream read_file("ActandAroSOC_SingleParticlePic_DataInput.cfg");
    assert(read_file.is_open());
    
    //Read in constants
    read_file >> Nsite;
    read_file >> Npart;
    read_file >> hoprat;
    read_file >> gamma;
    read_file >> phi_rat;
    read_file >> output;
    read_file >> otherout;
    read_file >> currentout;
    
    double Phi = phi_rat*Pi;
    //cout << "phi: " << Phi << endl;
    t_2 = t_1/hoprat;
    
    ofstream fout(output);
    assert(fout.is_open());
    fout.setf(ios::scientific);
    fout.precision(11);
    
    ofstream Fout(otherout);
    assert(Fout.is_open());
    Fout.setf(ios::scientific);
    Fout.precision(11);
    
    ofstream CurOut(currentout);
    assert(CurOut.is_open());
    CurOut.setf(ios::scientific);
    CurOut.precision(11);
    
    //Construct class
    Hamiltonian ham(t_1, t_2, gamma, Nsite);
    SlaterDet SD(ham);
    
    //Execute Functions
    //ham.Total_Ham(1, Phi);
    //ham.Peierls_Hamiltonian_pb();
    //ham.Diagonalize(fout, Fout);
    
    //Slater Determinant Functions
//    SD.SetNumber(Npart);
//    SD.FockBasis(ham);
//    SD.Get_Gstate_FB(ham);
    
    //Current
    //ProbCurrent(ham);
    
//    complex<double> J34_up = ham.BondCurrent(1, 2, 3);
//    complex<double> J34_dn = ham.BondCurrent(0, 8, 9);
//    cout << "Jup: " << J34_up << endl;
//    cout << "Jdn: " << J34_dn << endl;
//    
//    cout << "Js: " << J34_up - J34_dn << endl;
//    cout << "Jq: " << J34_up + J34_dn << endl;
    
//    for(int i = 0; i < 6; i++)
//    {
//        complex<double> J34_up = ham.BondCurrent(i, 2, 3);
//        
//        CurOut << i << " " << J34_up << endl;
//    }
    
    //CurrWFlux(ham, CurOut);
    
    //Eval over changing flux
    //EnergySpec(ham, fout);
    
    //Eval over gamma loop with current
    for(int g = 0; g <= Gm ; g++)
    {
        gamma = g*dg;
        cout << "g: " << g <<endl;
        ham.ResetGamma(gamma);
        ham.Total_Ham(1, Phi);
        //ham.Diagonalize(fout, Fout);
        
        VectorXd E = ham.ReturnEval();
        //complex<double> J34_up = ham.BondCurrent(0, 2, 3);
        for(int i = 0; i < 6; i++)
        {
            complex<double> Current = ham.BondCurrent(i, 2, 3);
            cout << gamma << " " << Current << endl;
            CurOut << gamma << " " << i << " " << Current.real() << endl;
        }
//        complex<double> J34_dn = ham.BondCurrent(1, 8, 9);
//        CurOut << gamma << " " << J34_up.real() - J34_dn.real() << " " << J34_up.real() + J34_dn.real() << endl;
        fout << gamma << " " << E.transpose() << endl;
        ham.ResetHam();
    }
    
    CurOut.close();
    fout.close();
    Fout.close();
    return 0;
}

void ProbCurrent(Hamiltonian &h)
{
    complex<double> J34_up;
    complex<double> J34_dn;
    complex<double> J_Tot;
    complex<double> JsT;
    
    J34_up = h.BondCurrent(0, 1, 2);
    J34_dn = h.BondCurrent(1, 7, 8);
    J_Tot = J34_dn + J34_up;
    JsT = J34_up - J34_dn;
    
    cout << "Jst: " << JsT << endl;
    cout << "Jup: " << J34_up << " Jdn: " << J34_dn << " Jqt: " << J_Tot << endl;
}

void EnergySpec(Hamiltonian &h, ofstream &out)
{
    double Phi_max = 1.0;
    double dp = .001;
    double P_it = Phi_max/dp;
    double phi;
    
    for(int p = 0; p <= P_it; p++)
    {
        //cout << p<< endl;
        phi = dp*p*Pi;
        cout << "Phi: " << phi << endl;
        h.Total_Ham(1, phi);
        VectorXd E = h.ReturnEval();
        //out << phi/Pi << " " << E(0) << " " << E(1) << endl;
        out << phi/Pi << " " << E.transpose() << endl;
    }
}

void CurrWFlux(Hamiltonian &h, ofstream &out)
{
    double Phi_max = 1.0;
    double dp = .01;
    double P_it = Phi_max/dp;
    double phi;
    
    complex<double> J34_up;
    complex<double> J34_dn;
    complex<double> J;
    double Jq;
    double Js;
    
    for(int p = 0; p <= P_it; p++)
    {
        //cout << p<< endl;
        phi = dp*p*Pi;
        cout << "Phi: " << phi << endl;
        h.Total_Ham(1, phi);
//        J34_up = h.BondCurrent(1, 2, 3);
//        J34_dn = h.BondCurrent(0, 8, 9);
//        Jq = J34_dn.real() + J34_up.real();
//        Js = J34_up.real() - J34_dn.real();
//        
//        cout << phi << " " << Js << " " << Jq << endl;
//        out << phi/Pi << " " << Js << " " << Jq << endl;
        for(int i = 0; i < 6; i++)
        {
        J = h.BondCurrent(i, 2, 3);
        cout << phi << " " << J << endl;
        out << phi/Pi << " " << i << " " << J.real() << endl;
        }
        h.ResetHam();

    }
 
}
