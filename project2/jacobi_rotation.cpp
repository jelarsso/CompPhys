#include "jacobi_rotation.hpp"
#include <cmath>
#include <iostream>

double off_diagonl(double** array, int size){
    double sum = 0;
    for (int i=0;i<size;i++){
        for (int j=0;j<size;j++){
            if (j==i){
                continue;
            }
            sum = sum + array[i][j]*array[i][j];
        }
    }
    return std::sqrt(sum);
};

void rotate(double** array, int size, int l, int k, double c, double s){
    double akk,all,akl,aik,ail;

    akk = array[k][k];
    all = array[l][l];
    akl = array[k][l];

    array[k][k] = akk*c*c  - 2*akl*c*s + all*s*s;
    array[l][l] = all*c*c + 2*akl*c*s + akk*s*s;
    array[k][l] = 0;
    array[l][k] = 0;

    for (int i=0;i<size;i++){
        if (i==k or i==l){
            continue;
        }
        aik = array[i][k];
        ail = array[i][l];

        array[i][k] = aik*c - ail*s;
        array[k][i] = array[i][k];
        array[i][l] = ail*c - aik*s;
        array[l][i] = array[i][l];
    }
};




int jacobi_rotate(double** array, double* eigvals, int size, double tolerance){
    /*

    Inputs:
    double** array: a 2d array of values of size nxn
    double* eigavls: an array of length n to store eigenvalues in
    int size: size n of the array
    double tolerance: the tolerance factor in the jacobi jacobi_rotation

    Output:
    int: the number of iterations that was needed
    */
    int number_of_iterations = 0;
    int k = 0;
    int l = 0;
    double tau;
    double t,c,s;

    while (off_diagonl(array,size)>tolerance){
        //find largest value:
        for (int i=0;i<size;i++){
            for (int j=0;j<size;j++){
                if (j==i){
                    continue;
                }
                if (std::abs(array[i][j])>std::abs(array[k][l])){
                    k = i;
                    l = j;
                }
            }
        }
        if (array[k][l]!=0){
        //find tau, t, c and s
        tau = (array[l][l] - array[k][k])/(2*array[k][l]);
        if (tau>=0){
            t = 1*(tau + std::sqrt(1 + tau*tau)); // or minus or something else, not sure!
        }else{
            t = -1*(-tau + std::sqrt(1 + tau*tau)); // or minus or something else, not sure!
        }
        c = 1/std::sqrt(1+t*t);
        s = t*c;
        }else{
            c = 1;
            s = 0;
        }

        //rotate
        rotate(array,size,l,k,c,s);
        number_of_iterations++;
    }
    for (int i=0;i<size;i++){
    eigvals[i] = array[i][i];
    }
return number_of_iterations;
};
