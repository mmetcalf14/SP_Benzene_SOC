//
//  main.cpp
//  ActandAro_SOC_SPM
//
//  Created by mekena McGrew on 1/9/17.
//  Copyright © 2017 Mekena Metcalf. All rights reserved.
//

#include <iostream>
#include "Hamiltonian.h"

void ProbCurrent(Hamiltonian &h);

int main(int argc, const char * argv[]) {
    
    
    //declare constants
    double t_1 = 1.0;
    double t_2;
    double hoprat;
    double gamma;
    double dg = .01;
    int Gm = 2/dg;
    
    int Nsite;
    int Npart;
    char output[60];
    char otherout[60];
    //call config file eventually
    ifstream read_file("ActandAroSOC_SingleParticlePic_DataInput.cfg");
    assert(read_file.is_open());
    
    //Read in constants
    read_file >> Nsite;
    read_file >> Npart;
    read_file >> hoprat;
    read_file >> gamma;
    read_file >> output;
    read_file >> otherout;
    
    t_2 = t_1/hoprat;
    
    ofstream fout(output);
    assert(fout.is_open());
    fout.setf(ios::scientific);
    fout.precision(11);
    
    ofstream Fout(otherout);
    assert(Fout.is_open());
    Fout.setf(ios::scientific);
    Fout.precision(11);
    
    //Construct class
    Hamiltonian ham(t_1, t_2, gamma, Nsite);
    
    //Execute Functions
    ham.Total_Ham();
    //ham.Peierls_Hamiltonian_pb();
    //ham.Diagonalize(fout, Fout);
    
    //Current
    //ProbCurrent(ham);
    
    complex<double> J34_up = ham.BondCurrent(0, 2, 3);
    complex<double> J34_dn = ham.BondCurrent(1, 8, 9);
    cout << "Js: " << J34_up - J34_dn << endl;
    
    //Eval over gamma loop
//    for(int g = 0; g <= Gm ; g++)
//    {
//        gamma = g*dg;
//        cout << "g: " << g <<endl;
//        ham.ResetGamma(gamma);
//        ham.Total_Ham();
//        ham.Diagonalize(fout, Fout);
//        ham.ResetHam();
//    }
    
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
