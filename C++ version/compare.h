//
//  compare.h
//  Polynomial Homotopy
//
//  Created by Dimitri Leggas on 10/22/14.
//
//

#ifndef __Polynomial_Homotopy__compare__
#define __Polynomial_Homotopy__compare__

#include <iostream>
#include <vector>
#include <complex>

using namespace std;

void sort(vector<complex<double>>*, int);
void conjugate(vector<complex<double>>*, int);
bool same(vector<complex<double>>, vector<complex<double>>, int);
void shift_origin(vector<complex<double>>, vector<complex<double>>*, int, int);
bool same_set(vector<complex<double>>, vector<complex<double>>, int);

#endif /* defined(__Polynomial_Homotopy__compare__) */
