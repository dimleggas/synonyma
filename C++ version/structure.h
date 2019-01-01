//
//  structure.h
//  Polynomial Homotopy v5
//
//  Created by Dimitri Leggas on 5/26/15.
//  Copyright (c) 2015 Dimitri Leggas. All rights reserved.
//

#ifndef __Polynomial_Homotopy_v5__structure__
#define __Polynomial_Homotopy_v5__structure__

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <complex>
#include "homotopy.h"
#include "intensity.h"

using namespace std;

class structure : public homotopy {
    
protected:
    const int atoms;
    int opt = 0;
    vector<complex<double>> z_final; //used to generate the system
    vector<complex<double>> I_final;
    vector<complex<double>> z_start; //generates any potential starting system
    vector<complex<double>> I_start;
    intensity int_calc;
    
    void hh(vector<complex<double>>*, vector<complex<double>>, double, int);
    void jj(vector<complex<double>>*, vector<complex<double>>, double, int);
    
    vector<complex<double>> old_solution;
    void store_solution();
    
    /* Solving */
    void newton_corrector(int, double); // overriding from homotopy
    void naive(int); // overriding from homotopy
    void track();
    void random_restart(); // randomly resets starting system
    void correlated_intensity_restart(); // resets starting system so that intensities  positively correlate with target
    void coordinate_reset(); // randomly resets one coordinate
    bool check_solution();
    bool check_intensity(); // should not use
    
    
    /* Analysis */
    int attempts = 0; // number or restarts until a solution is found
    int meaning = 0; // number or restarts until a solution is found
    bool swap = false; // did the solution switch between meaningful/non-meaningful
    
    int order_violation(vector<double>, vector<double>, int);
    void analyze_intensities();
    
public:
    structure(int num, double time_step): homotopy(3 * (num - 1), time_step), atoms(num), int_calc(intensity(num)) {
        int_calc.generate_coordinates(&z_final);
        int_calc.generate_coordinates(&z_start);
        int_calc.generate_intensities(z_final, &I_final);
        int_calc.generate_intensities(z_start, &I_start);
        
        for(int i = 0; i < n; i++)
            z_current[i] = z_start[i];
        
        old_solution.resize(5 * n);
    };
    
    void generate_intensities();
    void set_intensities(vector<complex<double>>);
    void set_z(vector<complex<double>>, int opt = 0);
    void restart(int);
    void solve(int);
    void recursive_solve(int);
    void print_state();
    bool meaningful();
    
    /* set/get */
    int get_attempts();
    int get_meaningful();
    
};


#endif /* defined(__Polynomial_Homotopy_v5__structure__) */
