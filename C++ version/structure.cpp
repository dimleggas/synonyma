//
//  structure.cpp
//  Polynomial Homotopy v5
//
//  Created by Dimitri Leggas on 5/26/15.
//  Copyright (c) 2015 Dimitri Leggas. All rights reserved.
//

#include <ctime>
#include "structure.h"
#include "homotopy.h"
#include "matrix.h"
#include "compare.h"
#include "data.h"

/*
 * Homotopy methods
 */

// homotopy
void structure::hh(vector<complex<double>>* homo, vector<complex<double>> z, double t, int dim){
    complex<double> t_final (t, 0.0);
    complex<double> t_start (1 - t, 0.0);
    
    for(int i = 0; i < dim; i++)
        (*homo)[i] = int_calc.system(0, z, i, 0, opt, t) - t_start * I_start[i] - t_final * I_final[i];
}

// jacobian
void structure::jj(vector<complex<double>>* jaco, vector<complex<double>> z, double t, int dim){
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            (*jaco)[element(i, j, dim)] = int_calc.system(1, z, i, j, opt, t);
}

/*
 * Override the newton_corrector to call structure specific homotopy and jacobian functions
 * Ovveride the naive to utilize restart
 */

// corrects the homotopy path
void structure::newton_corrector(int max_iteration, double cutoff) {
    int count =  0;
    bool cont; // continue?
    vector<complex<double>> z_last (n);
    
    do{
        cont = false;
        
        // calculate current homotopy and its Jacobian
        hh(&homo, z_current, t, n);
        jj(&jacobian, z_current, t, n);
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

// follow the homotopy path (naive method => constant dt)
void structure::naive(int type) {
    
    track();
    
    //analyze_intensities();
    if (meaningful())
        meaning++;

    
    if (!meaningful() || !check_solution() || !same_set(z_current, z_final, atoms)) {
        restart(type);
        naive(type);
    }
}

// follow the homotopy path (naive method => constant dt)
void structure::track() {
    int count = 0;
    int it = 1 / dt;
    
    t = 0;
    while (count < it) {
        count++;
        t+=dt;
        structure::newton_corrector(15, 10E-6);
        
        for(int i = 0; i < n; i++)
            if(abs(z_current[i]) > 100) return;
        
        // see if swapped or if need to reset
        if (!meaningful())
            swap = true;
    }
}

/* Setting up the system */
void structure::generate_intensities(){
    int_calc.generate_coordinates(&z_final);
    int_calc.generate_intensities(z_final, &I_final);
}

void structure::set_intensities(vector<complex<double>> i_list){
    for(int i = 0; i < n; i++)
        I_final[i] = i_list[i];
}

void structure::set_z(vector<complex<double>> z, int opt){
    for (int i = 0; i < n; i++){
        z_current[i] = z[i];
        z_start[i] = z[i];
    }
    int_calc.generate_intensities(z_start, &I_start);
}

/* 
 * For solve() and restart()
 * 0 for normal restart
 * 1 for correlated intensity pick
 * 2 for matching lower order intensities 
 * 3 ensures both 1 and 2
 */

void structure::solve(int type = 0){
    attempts = 0;
    meaning = 0;

    restart(type);
    naive(type);
}

void structure::restart(int type){
    swap = false;
    attempts++;
    
    // random restart
    if (0 == type) {
        int_calc.generate_coordinates(&z_start);
        int_calc.generate_intensities(z_start, &I_start);
    }
    
    // correlated intensitity restart
    if (1 == type) {
        bool smart = false;
        do{
            int_calc.generate_coordinates(&z_start);
            int_calc.generate_intensities(z_start, &I_start);
            
            if (real(correlation(I_start, I_final, n)) > 0.5)
                smart = true;
        } while (!smart);
    }
    
    // matching lower order intensities
    if (2 == type) {
        complex<double> one (1, 0);
        complex<double> two (2, 0);
        complex<double> four (4, 0);
        int_calc.generate_coordinates(&z_start);
        // change three coordinates so that lowest 3 intensities match final
        for (int i = 0; i < 3; i++){
            complex<double> a (1, 0);
            complex<double> b (2 - real(I_final[i]), 0);
            complex<double> c (1, 0);
            for (int j = 1; j < atoms - 1; j++){
                a += one / z_start[element(j, i, 3)];
                b += z_start[element(j, i, 3)] + one / z_start[element(j, i, 3)];
                for (int k = 1; k < atoms - 1; k++)
                    b += z_start[element(j, i, 3)] / z_start[element(k, i, 3)];
                c += z_start[element(j, i, 3)];
            }
            z_start[i] = (-b + sqrt(pow(b, 2) - four * a * c)) / (two * a);
        }
        int_calc.generate_intensities(z_start, &I_start);
    }
    
    // matching lower order intensities and postively correlated intensities
    if (3 == type) {
        bool smart = false;
        do{
            complex<double> one (1, 0);
            complex<double> two (2, 0);
            complex<double> four (4, 0);
            int_calc.generate_coordinates(&z_start);
            // change three coordinates so that lowest 3 intensities match final
            for (int i = 0; i < 3; i++){
                complex<double> a (1, 0);
                complex<double> b (2 - real(I_final[i]), 0);
                complex<double> c (1, 0);
                for (int j = 1; j < atoms - 1; j++){
                    a += one / z_start[element(j, i, 3)];
                    b += z_start[element(j, i, 3)] + one / z_start[element(j, i, 3)];
                    for (int k = 1; k < atoms - 1; k++)
                        b += z_start[element(j, i, 3)] / z_start[element(k, i, 3)];
                    c += z_start[element(j, i, 3)];
                }
                z_start[i] = (-b + sqrt(pow(b, 2) - four * a * c)) / (two * a);
            }
            int_calc.generate_intensities(z_start, &I_start);
            
            if (real(correlation(I_start, I_final, n)) > 0.5)
                smart = true;
        } while (!smart);
    }
    
    for(int i = 0; i < n; i++)
        z_current[i] = z_start[i];
    
}

void structure::recursive_solve(int op){
    attempts = 0;
    meaning = 0;
    opt = op;
    
    vector<complex<double>> z_sub (6);
    vector<complex<double>> z_copy (n);
    vector<complex<double>> I_sub (6);

    for(int at = 3; at < atoms; at++){
        for(int i = 0; i < z_sub.size(); i++)
            z_copy[i] = z_sub[i];
        
        z_sub.resize(3 * (at - 1));
        I_sub.resize(3 * (at - 1));
        for (int i = 0; i < z_sub.size(); i++)
            z_copy[i] = z_sub[i];
        
        complex<double> div (pow(n / at, 2), 0);
        for (int i = 0; i < 3 * (at - 1); i++)
            I_sub[i] = I_final[i] / div;
        
        structure subsolver(at, 1E-3);
        subsolver.set_intensities(I_sub);
        
        for (int i = 0; i < z_sub.size(); i++)
            z_copy[i] = z_sub[i];
        
        // decide whether to set initial guess or augment guess
        int begin = 4 * (at - 2);
        if(at == 4) begin = 0;
        int count = 0;
        do{
            for (int i = 0; i < 3 * (at - 1); i++){
                if (i < begin)
                    z_sub[i] = z_copy[i];
                else{
                    double x = 2 * M_PI * int_calc.get_random();
                    complex<double> temp (cos(x), sin(x));
                    z_sub[i] = temp;
                }
            }
            subsolver.set_z(z_sub, 1);
            subsolver.track();
            //subsolver.print_state();
            count++;
            if (count > 20){
                cout << "{0,0},";
                return;
            }
            subsolver.get_z(&z_sub);
        }while(!subsolver.meaningful());
        subsolver.get_z(&z_sub);
    }
    
    do{
        // augment the solution for N-1 to solve the system for N atoms
        for (int i = 0; i < n - 3; i++)
            z_start[i] = z_sub[i];
        for (int i = n - 3; i < n; i++){
            double x = 2 * M_PI * int_calc.get_random();
            complex<double> temp (cos(x), sin(x));
            z_start[i] = temp;
        }
        int_calc.generate_intensities(z_start, &I_start, opt, 0);
        cout << "Target/Start intensities" << endl;
        for (int i = 0; i < n; i++)
            cout << I_final[i] << " " << I_start[i] << endl;
        
        for(int i = 0; i < n; i++)
            z_current[i] = z_start[i];
        
        cout << "solving target system " << endl;
        for(int i = 0; i < n; i++)
            cout<< z_start[i] << endl;
        
        track();
        //print_state();
        attempts++;
        if (meaningful())
            meaning++;
    }while(!same_set(z_current, z_final, atoms));
    cout << "{" << attempts << "," << meaning << "},";
    hh(&homo, z_current, t, n);
    for(int i = 0; i < n; i++)
        cout << homo[i] << endl;
    
}

/*
 * Set/Get/Print
 */

int structure::get_attempts(){
    return attempts;
}

int structure::get_meaningful(){
    return meaning;
}

void structure::print_state(){
    if(check_solution()){
        cout << "Solution found";
        if (meaningful()) cout << " and is meaningful." << endl;
        else cout << ", but is not meaningful." << endl;
        cout << "It took " << attempts << " homotopy tracks." << endl;
    }
    else{
        if (meaningful()) cout << "Current homotopy solution is meaningful." << endl;
        else cout << "Current homotopy solution is not meaningful" << endl;
    }
}

/*
 * Diagnostics and Analysis
 */

// SOLUTION/SYSTEM DIAGNOSTICS

// determine if current solution is solution to target system
bool structure::check_solution() {
    for(int i = 0; i < n; i++)
        if(abs(int_calc.system(0, z_current, i, 0, 0, 0) - I_final[i]) > 1E-6) return false;
    return true;
}

// determine if current solution is meaningful
bool structure::meaningful() {
    for(int i = 0; i < n; i++)
        if(abs(abs(z_current[i]) - 1) > 1E-6) return false;
    return true;
}

// a check to see if two solutions belong to the same set -- use the more rigorous method in compare
bool structure::check_intensity() {
    for(int i = n; i < n + 3; i++)
        if(abs(int_calc.system(0, z_final, i, 0, 0, 0) - int_calc.system(0, z_current, i, 0, 0, 0)) > 1E-6) return false;
    return true;
}

// track previous solutions
void structure::store_solution(){
    for (int i = 4; i > 0; i--)
        for (int j = 0; j < n; j++)
            old_solution[element(i, j, n)] = old_solution[element(i - 1, j, n)];
    for (int i = 0; i < n; i++)
        old_solution[element(0, i, n)] = z_current[i];
}

// INTENSITY ANALYSIS

// calculate the number of order violations between two intensity vectors
int structure::order_violation(vector<double> v, vector<double> w, int n){
    int count = 0;
    for (int i = 0; i < n - 1; i++)
        for (int j = i + 1; j < n; j++)
            if ((v[i] < v[j]) != (w[i] < w[j])) count++;
    return count;
}

// analyze relationship between starting and ending intensities
void structure::analyze_intensities() {
    if(!meaningful())
        cout << "{0,";
    else {
        if (swap)
            cout << "{1,";
        else
            cout << "{2,";
    }
    vector<double> Is (n);
    vector<double> If (n);
    for (int i = 0; i < n; i++){
        Is[i] = real(I_start[i]);
        If[i] = real(I_final[i]);
    }
    cout << order_violation(Is, If, n) << ",";
    cout << real(correlation(I_start, I_final, n)) << "},"; // use real part since imag ~ 1E-17 (roundoff)
}


/*
 * Old methods
 */

void structure::random_restart(){
    swap = false;
    attempts++;
    int_calc.generate_coordinates(&z_start);
    int_calc.generate_intensities(z_start, &I_start);
    for(int i = 0; i < n; i++)
        z_current[i] = z_start[i];
}

void structure::correlated_intensity_restart(){
    swap = false;
    attempts++;
    
    bool smart = false;
    do{
        int_calc.generate_coordinates(&z_start);
        int_calc.generate_intensities(z_start, &I_start);
        
        if (real(correlation(I_start, I_final, n)) > 0.5)
            smart = true;
    } while (!smart);
    
    for(int i = 0; i < n; i++)
        z_current[i] = z_start[i];
}

void structure::coordinate_reset(){
    swap = false;
    
    for(int i = 0; i < n; i++)
        z_current[i] = old_solution[element(4, i, n)];
    
    double r = int_calc.get_random();
    complex<double> temp (cos(2 * M_PI * r), sin(2 * M_PI * r));
    z_current[0] = temp;
    
    int_calc.generate_intensities(z_current, &I_start);
}

