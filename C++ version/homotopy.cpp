//
//  homotopy.cpp
//  Polynomial Homotopy v5
//
//  Created by Dimitri Leggas on 6/03/15.
//  Copyright (c) 2015 Dimitri Leggas. All rights reserved.
//

#include <iostream>
#include "homotopy.h"
#include "matrix.h"

using namespace std;

/* Set/Get */
void homotopy::set_system(vector<complex<double>> z_start, compute seth, compute setj){
    for (int i = 0; i < n; i++)
        z_current[i] = z_start[i];
    set_homotopy = seth;
    set_jacobian = setj;
}

void homotopy::set_time(double time) {
    if (time < 0 || time > 1) return;
    dt = time;
}

void homotopy::get_z(vector<complex<double>>* z){
    for (int i = 0; i < n; i++)
        (*z)[i] = z_current[i];
}

/* Homotopy tracking */

// estimate next solution on homotopy path
void homotopy::linear_step(){}

// corrects the homotopy path
void homotopy::newton_corrector(int max_iteration, double cutoff) {
    int count =  0;
    bool cont; // continue?
    vector<complex<double>> z_last (n);
    
    do{
        cont = false;
        
        // calculate current homotopy and its Jacobian
        set_homotopy(&homo, z_current, t, n);
        set_jacobian(&jacobian, z_current, t, n);
        inverse(jacobian, &jacobian_inverse, n);

        // properly adjust the solution
        matrix_product(jacobian_inverse, n, n, homo, n, 1, &dz);
        for (int i = 0; i < n; i++)
            z_current[i] -= dz[i];
        
        // compare previous and new solution
        for (int i = 0; i < n; i++) {
            if (abs(dz[i]) > cutoff){
                cont = true;
                break;
            }
        }

        
        count++;
    }while(cont && count < max_iteration);
    
}

/* Solve */

// take one iteration down the homotopy path
// random solve by default
void homotopy::naive() {
    int count = 0;
    int it = 1 / dt + 1;
    
    t = -dt;
    while (count < it) {
        count++;
        t+=dt;
        
        // ensure not diverging
        for (int i = 0; i < n; i++)
            if (abs(z_current[i]) > 1000) return;
    
        //linear_step();
        newton_corrector(15, 10E-6);
    }
}

void homotopy::step() {
    t+=dt;
    
    // ensure not diverging
    for (int i = 0; i < n; i++)
        if (abs(z_current[i]) > 1000) return;
    
    newton_corrector(15, 10E-6);
}

/* Testing */
void homotopy::test(){
    
}

