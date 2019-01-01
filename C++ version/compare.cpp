//
//  compare.cpp
//  Polynomial Homotopy
//
//  Created by Dimitri Leggas on 5/27/15.
//
//

#include <complex>
#include "compare.h"
#include "matrix.h"

using namespace std;

bool same_set(vector<complex<double>> z_current, vector<complex<double>> z_final, int atoms) {
    int n = 3 * (atoms - 1);
    
    // sort the coordinates used to generate the system and the found solution
    sort(&z_current, atoms);
    sort(&z_final, atoms);
    
    // copy the solution to "shifted" array
    vector<complex<double>> shifted (n);
    for (int i = 0; i < n; i++)
        shifted[i] = z_current[i];
    
    // compare solution and conjugate solution
    if (same(z_final, shifted, n))
        return true;
    conjugate(&shifted, n);
    sort(&shifted, atoms);
    if (same(z_final, shifted, n))
        return true;
    
    // if necessary, shift origin to each atom
    for (int i = 0; i < atoms - 1; i++) {
        shift_origin(z_current, &shifted, atoms, i);
        sort(&shifted, atoms);
        if (same(z_final, shifted, n))
            return true;
        conjugate(&shifted, n);
        sort(&shifted, atoms);
        if (same(z_final, shifted, n))
            return true;
    }
    
    // if no match, return false
    return false;
}

// sort first based on real portion of first coordinate, then imagninary
void sort(vector<complex<double>>* coords, int atoms) {
    bool swap = false;
    do {
        for (int i = 0; i < atoms - 2; i++) {
            
            // decide whether or not to swap
            swap = false;
            if (abs(real((*coords)[element(i, 0, 3)]) - real((*coords)[element(i + 1, 0, 3)])) > 1E-6) {
                if (real((*coords)[element(i, 0, 3)]) - real((*coords)[element(i + 1, 0, 3)]) > 0.)
                    swap = true;
            }
            else {
                if (imag((*coords)[element(i, 0, 3)]) - imag((*coords)[element(i + 1, 0, 3)]) > 0.)
                    swap = true;
            }
            
            // swap
            if (swap){
                for (int j = 0; j < 3; j++) {
                    complex<double> temp = (*coords)[element(i, j, 3)];
                    (*coords)[element(i, j, 3)] = (*coords)[element(i + 1, j, 3)];
                    (*coords)[element(i + 1, j, 3)] = temp;
                }
            }
        }
    }while (swap);
}

void shift_origin(vector<complex<double>> coords, vector<complex<double>>* shifted, int atoms, int origin) {
    vector<complex<double>> shift (3);
    for (int i = 0; i < 3; i++)
        shift[i] = conj(coords[element(origin, i, 3)]);
    for (int i = 0; i < atoms - 1; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == origin)
                (*shifted)[element(i, j, 3)] = shift[j];
            else
                (*shifted)[element(i, j, 3)] = coords[element(i, j, 3)] * shift[j];
        }
    }
}

// compare equality of matrices within error of 1E-6
bool same(vector<complex<double>> A, vector<complex<double>> B, int n) {
    for (int i = 0; i < n; i++)
        if(abs(A[i]-B[i]) > 1E-6) return false;
    return true;
}

// take the conjugate of every element in the matrix
void conjugate(vector<complex<double>>* A, int n) {
    for (int i = 0; i < n; i++)
        (*A)[i] = conj((*A)[i]);
}

