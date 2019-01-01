//
//  data.h
//  Polynomial Homotopy v5
//
//  Created by Dimitri Leggas on 5/28/15.
//  Copyright (c) 2015 Dimitri Leggas. All rights reserved.
//

#ifndef __Polynomial_Homotopy_v5__data__
#define __Polynomial_Homotopy_v5__data__

#include <stdio.h>
#include <vector>
#include <complex>

using namespace std;

complex<double> mean(vector<complex<double>>, int );
complex<double> standard_deviation(vector<complex<double>>, int);
complex<double> covariance(vector<complex<double>>, vector<complex<double>>, int);
complex<double> correlation(vector<complex<double>>, vector<complex<double>>, int);


#endif /* defined(__Polynomial_Homotopy_v5__data__) */