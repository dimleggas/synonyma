//
//  data.cpp
//  Polynomial Homotopy v5
//
//  Created by Dimitri Leggas on 5/28/15.
//  Copyright (c) 2015 Dimitri Leggas. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <vector>
#include <complex>
#include "data.h"
#include "matrix.h"

using namespace std;

complex<double> mean(vector<complex<double>> v, int n){
    complex<double> sum (0, 0);
    for (int i = 0; i < n; i++)
        sum += v[i];
    complex<double> div (n, 0);
    return sum / div;
}

complex<double> standard_deviation(vector<complex<double>> v, int n){
    complex<double> m =  mean(v, n);
    complex<double> sum (0, 0);
    for (int i = 0; i < n; i++)
        sum += pow(v[i] - m, 2);
    complex<double> div (n - 1, 0);
    sum /= div;
    return sqrt(sum);
}

complex<double> covariance(vector<complex<double>> v, vector<complex<double>> w, int n){
    complex<double> mv = mean(v, n);
    complex<double> mw = mean(w, n);
    
    complex<double> sum (0, 0);
    for (int i = 0; i < n; i++)
        sum += (v[i] - mv) * (w[i] - mw);
    
    complex<double> div (n - 1, 0);
    return sum /div;
}

complex<double> correlation(vector<complex<double>> v, vector<complex<double>> w, int n){
    complex<double> sdv = standard_deviation(v, n);
    complex<double> sdw = standard_deviation(w, n);

    return covariance(v, w, n) / (sdv * sdw);
}

