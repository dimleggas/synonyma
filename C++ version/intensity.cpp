//
//  intensity.cpp
//  Polynomial Homotopy v5
//
//  Created by Dimitri Leggas on 5/26/15.
//  Copyright (c) 2015 Dimitri Leggas. All rights reserved.
//

#include <iostream>
#include <complex>
#include "intensity.h"

using namespace std;


/*
    process = 0 for equation calculation
    otherwise derivative calculation
    derivative only used if process != 0
 
    opt !=0 for perturbed system
*/
complex<double> intensity::system(int process, vector<complex<double>> z, int eqn, int derivative, int opt, double t){
    complex<double> zero (0.0, 0.0);
    int *intensity_list = hkl[eqn];
    
    complex<double> one (1, 0);
    complex<double> nn (n, 0);
    complex<double> tt (t, 0);
    
    // set contribution of origin (depending on if perturbed or not)
    complex<double> F (1, 0.0);
    complex<double> Fconj (1, 0.0);
    // reassign if necessary
    if (0 != opt){
        F = (nn - tt) / (nn - one);
        Fconj = (nn - tt) / (nn - one);
    }
    // sum up contribution of atoms to the structure factor
    for (int i = 0; i < atoms - 1; i++) {
        complex<double> temp1 (1.0, 0.0);
        complex<double> temp2 (1.0, 0.0);
        // if perturbed system, weight the atoms properly
        if (0 != opt){
            if (i != atoms - 2) {
                temp1 *= (nn - tt) / (nn - one);
                temp2 *= (nn - tt) / (nn - one);
            }
            else{
                temp1 *= tt;
                temp2 *= tt;
            }
        }
        // calculate component contribution in the product
        for (int j = 0; j < 3; j++) {
            temp1 *= pow(z[3 * i + j], intensity_list[j]);
            temp2 *= pow(z[3 * i + j], -intensity_list[j]);
        }
        
        F += temp1;
        Fconj += temp2;
    }
    
    if (process == 0) return F * Fconj;
    
    // find the term where the differentiating variable is located, which factor it is
    int term = derivative / 3;
    int fact = derivative % 3;
    
    // if the term is constant return 0
    if(intensity_list[fact] == 0) {
        return zero;
    }
    else {
        complex<double> dF (1.0, 0.0);
        complex<double> dFconj (1.0, 0.0);
        // if system is a perturbed subsystem
        if (0 != opt){
            if (term != atoms - 2) {
                dF *= (nn - tt) / (nn - one);
                dFconj *= (nn - tt) / (nn - one);
            }
            else{
                dF *= tt;
                dFconj *= tt;
            }
        }
        for (int i = 3 * term; i < 3 * (term + 1); i++) {
            if (i % 3 == fact) {
                complex<double> mult1 (intensity_list[i % 3], 0.0);
                complex<double> mult2 (-intensity_list[i % 3], 0.0);
                dF *= mult1 * pow(z[i], intensity_list[i % 3] - 1);
                dFconj *= mult2 * pow(z[i], -intensity_list[i % 3] - 1);
            }
            else {
                dF *= pow(z[i], intensity_list[i % 3]);
                dFconj *= pow(z[i], -intensity_list[i % 3]);
            }
        }
        return dF * Fconj + F * dFconj;
    }
}

// randomly generates coordinates
void intensity::generate_coordinates(vector<complex<double>>* z){
    (*z).resize(n);
    
    // generate a solution to the system
    for (int i = 0; i < n; i++) {
        r = (double) rand() / RAND_MAX;
        complex<double> temp (cos(2 * M_PI * r), sin(2 * M_PI * r));
        (*z)[i] = temp;
    }
}

// calculates intensities from randomly generated coordinates
void intensity::generate_intensities(vector<complex<double>> z, vector<complex<double>>* I, int opt, int t){
    (*I).resize(n);
    
    // calculate intensities for target system
    for (int i = 0; i < n; i++)
        (*I)[i] = system(0, z, i, 0, opt, t);
}


double intensity::get_random(){
    return(double) rand() / RAND_MAX;
}

