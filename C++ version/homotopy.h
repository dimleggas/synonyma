//
//  homotopy.h
//  Polynomial Homotopy v5
//
//  Created by Dimitri Leggas on 5/22/15.
//  Copyright (c) 2015 Dimitri Leggas. All rights reserved.
//

#ifndef __Polynomial_Homotopy_v5_1__homotopy__
#define __Polynomial_Homotopy_v5_1__homotopy__

#include <stdio.h>
#include <complex>
#include <vector>

using namespace std;

typedef void (*compute)(vector<complex<double>>*, vector<complex<double>>, double, int);

class homotopy{
    
protected:
    const int n;
    double t = 0;
    double dt;
    vector<complex<double>> z_current;
    vector<complex<double>> homo;
    vector<complex<double>> jacobian;
    vector<complex<double>> jacobian_inverse;
    vector<complex<double>> dz; // change in solution
    
    /* Functions */
    compute set_homotopy;  // sets the homotopy
    compute set_jacobian;  // sets the jacobian
    void linear_step();  // one homotopy step (1st order)
    void newton_corrector(int, double);   // corrects the homotopy estimate
    
public:
    homotopy(int dimension, double time_step): n(dimension), dt(time_step) {
        // intialize arrays
        z_current.resize(n);
        jacobian.resize(n * n);
        jacobian_inverse.resize(n * n);
        homo.resize(n);
        dz.resize(n);
    };
    
    void set_system(vector<complex<double>>, compute, compute);
    void set_time(double);
    void get_z(vector<complex<double>>*);
    void get_f(vector<complex<double>>*);
    void naive(); // naive homotopy
    void track(); // track homotopy
    void step(); //one homotopy step
    
    /* Clean-up */
    void clean_double(complex<double>**, int);
    void clean();
    
    /* Testing */
    void test();
};

#endif /* defined(__Polynomial_Homotopy_v5__homotopy__) */
