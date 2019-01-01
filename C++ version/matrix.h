//
//  matrix.h
//  Polynomial Homotopy
//
//  Created by Oleg Tsodikov on 6/3/14.
//  Modified by Dimitri Leggas on 5/26/15
//

#ifndef __Polynomial_Homotopy_v5__matrix__
#define __Polynomial_Homotopy_v5__matrix__

#include <complex>
#include <vector>

using namespace std;

int element(int, int, int);
void print_matrix(vector<complex<double>>, int);
void matrix_product(vector<complex<double>>, int, int, vector<complex<double>>, int, int, vector<complex<double>>*);
void inverse(vector<complex<double>>, vector<complex<double>>*, int);


#endif /* defined(__Polynomial_Homotopy_v5__matrix__) */
