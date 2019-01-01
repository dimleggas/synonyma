//
//  main.cpp
//  Polynomial Homotopy v5
//
//  Created by Dimitri Leggas on 6/03/15.
//  Copyright (c) 2015 Dimitri Leggas. All rights reserved.
//

#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include "structure.h"
#include "homotopy.h"
#include "matrix.h"
#include "data.h"
#include "intensity.h"


using namespace std;

int main(int argc, const char * argv[]) {
    
    vector<complex<double>> z (3);
    structure solver (2, 1E-3);
    cout << "{";
    for (int i = 0; i < 100; i++){
        solver.generate_intensities();
        solver.solve(0);
        cout << "{" << solver.get_attempts() << "," << solver.get_meaningful() << "},";
        solver.get_z(&z);
        for (int i = 0; i < 3; i++)
            cout << z[i] << " " << abs(z[i]) << endl;

//        solver.recursive_solve(1);
    }
    cout << "}";
//    solver.get_z(&z);
//    for (int i = 0; i < 6; i++)
//        cout << z[i] << " " << abs(z[i]) << endl;
    
    return 0;
}
