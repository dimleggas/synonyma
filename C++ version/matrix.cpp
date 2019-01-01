//
//  matrix.h
//  Polynomial Homotopy
//
//  Created by Oleg Tsodikov on 6/3/14.
//  Modified by Dimitri Leggas on 5/26/15
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <math.h>
#include <vector>

using namespace std;

int element(int row, int column, int row_length){
    return row * row_length + column;
}

void print_matrix(vector<complex<double>> v, int row_length){
    long num_rows = v.size() / row_length;
    for (int i = 0; i < num_rows; i++){
        for (int j = 0; j < row_length; j++)
            cout << v[element(i, j, row_length)];
        cout << endl;
    }
}

void matrix_product(vector<complex<double>> A, int s, int t, vector<complex<double>> B, int r, int q, vector<complex<double>>* C){
    if (t != r) return; // ensure matrices can be multiplied
    
    (*C).resize(s * q);
    
    complex<double> sum;
    complex<double> zero (0.0, 0.0);
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < q; j++) {
            sum = zero;
            for (int k = 0; k < t; k++)
                sum += A[element(i, k, t)] * B[element(k, j, q)];
            (*C)[element(i, j, q)] = sum;
        }
    }
}

void inverse(vector<complex<double>> matr, vector<complex<double>>* invmatr, int n){

    int i, j, pivot, pivotflag, pivoti = 0, nzero;
    complex<double> coeff, temp;
    complex<double> one(1.,0.), imagone(0.,1.);
    complex<double> zero(0.,0.);
    
    //initializing invmatr as an identity matrix
    (*invmatr).resize(n * n);
    for(i = 0; i < n; i++)
        for(j = 0; j < n; j++)
            if(j==i) (*invmatr)[element(i, j, n)] = one;
            else (*invmatr)[element(i, j, n)]=zero;
    
    //converting to zeros below the diagonal
    
    for(pivot = 0; pivot < n; pivot++){  // looking for the pivot row

        pivotflag = 0;
        
        for(i = 0; i < n; i++){
    
            { for(j = 0, nzero = 0; j < pivot; j++)
                if(matr[element(i, j, n)]==zero) nzero++;
                
                if(nzero==pivot && matr[element(i, pivot, n)]!=zero){
                    pivotflag = 1;  //found the non-zero element
                    pivoti = i;  //pivot row
                    break;
                }
            }
        }
        if(pivotflag == 0) { printf("\nSingular matrix, exiting...");
            exit(0);
        }
        
        //pivot row is found (pivoti), now swap this into the right place (pivot)
        
        if(pivoti != pivot)
            for(j = 0; j < n; j++){
                temp = matr[element(pivot, j, n)]; // save element
                matr[element(pivot, j, n)] = matr[element(pivoti, j, n)]; //replace element
                matr[element(pivoti, j, n)] = temp; //replace element
                
                // same for the augmented matrix
                temp = (*invmatr)[element(pivot, j, n)]; // save element
                (*invmatr)[element(pivot, j, n)] = (*invmatr)[element(pivoti, j, n)]; //replace element
                (*invmatr)[element(pivoti, j, n)] = temp; //replace element
            }
        
        
        //now substracting the pivot row multiplied by a correct scalar from the other rows to get zeros below diagonal
        
        
        for(i = pivot + 1; i < n; i++){
            if(matr[element(i, pivot, n)]!=zero){
                
                coeff = matr[element(i, pivot, n)] / matr[element(pivot, pivot, n)];
                
                for(j = 0; j < n; j++){
                    
                    if(matr[element(pivot, j, n)]!=zero)
                        matr[element(i, j, n)] -= matr[element(pivot, j, n)]*coeff;
                    
                    if(j==pivot)        /*elements below the pivot point turn into (0,0)*/
                        matr[element(i, j, n)]=zero;
                    
                    if((*invmatr)[element(pivot, j, n)]!=zero)
                        (*invmatr)[element(i, j, n)] -= (*invmatr)[element(pivot, j, n)]*coeff;
                }
            }
        }
        
    }
    
    // now all elements below diagonal should be 0
    
    // normalizing diagonal elements first, to simplify job
    
    for(i = 0; i < n; i++){
        for(j = i + 1; j < n; j++)
            matr[element(i, j, n)] /= matr[element(i, i, n)];
        
        for(j = 0; j < n; j++)  //need to do this for all elements of the row for invmatr
            (*invmatr)[element(i, j, n)] /= matr[element(i, i, n)];
        
        matr[element(i, i, n)] = one;  // turning the diagonal elements into 1 (not really needed)
    }
    
    
    // now turning elements above diagonal to 0, bottom to top
    for(pivot = n-1; pivot > 0; pivot--){
        for(i = pivot - 1; i >= 0; i--)
            if(matr[element(i, pivot, n)] != zero){
                coeff=matr[element(i, pivot, n)];
                for(j = 0; j < n; j++){
                    matr[element(i, j, n)] -= matr[element(pivot, j, n)] * coeff;
                    (*invmatr)[element(i, j, n)] -= (*invmatr)[element(pivot, j, n)] * coeff;
                }
            }
    }

}
