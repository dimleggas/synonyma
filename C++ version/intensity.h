//
//  intensity.h
//  Polynomial Homotopy v5
//
//  Created by Dimitri Leggas on 5/26/15.
//  Copyright (c) 2015 Dimitri Leggas. All rights reserved.
//

#ifndef __Polynomial_Homotopy_v5__intensity__
#define __Polynomial_Homotopy_v5__intensity__

#include <stdio.h>
#include <complex>
#include <vector>

using namespace std;

class intensity {
    const int atoms;    // number of atoms
    int n;    // system size
    double r; // random x in R intersection [0, 1]
    int hkl[39][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 1, 0 },
        { 1, 0, 1 }, { 0, 1, 1 }, { 1, 0, -1 }, { 1, -1, 0 }, { 0, 1, -1 },
        { 2, 1, 0 }, { 2, 0, 1 }, { 0, 2, 1 }, { 1, 2, 0 },
        { 0, 1, 2 }, { 1, 0, 2 }, {1, 3, 0} , {3, 1, 0} , {1, 0, 3}, {3, 0, 1},
        {0, 3, 1}, {0, 1, 3}, {3, 2, 0}, {3, 0, 2}, {2, 3, 0}, {2, 0, 3}, {0, 2, 3},
        {0, 3, 2}, {4, 1, 0}, {4, 0, 1}, {1, 4, 0}, {1, 0, 4}, {0, 1, 4}, {0, 4, 1},
        {5, 0, 1}, {5, 1, 0}, {5, 2, 0}, {5, 0, 2}, {1, 0, 5}, {1, 5, 0}
    };
    
    complex<double> term(int type);
    complex<double> dterm(int type);
    
public:
    intensity(int num): atoms(num) {
        n = 3 * (atoms - 1);
        
        srand (static_cast<unsigned int>(time(NULL)));
        r = (double) rand() / RAND_MAX;
    };

    complex<double> system(int, vector<complex<double>>, int, int, int, double); // calculates equation or its derivative
    void generate_coordinates(vector<complex<double>>*); // generate intensities
    void generate_intensities(vector<complex<double>>, vector<complex<double>>*, int opt = 0, int t = 0); // generate intensities
    
    double get_random();

};

#endif /* defined(__Polynomial_Homotopy_v5__intensity__) */
