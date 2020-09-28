#include "jacobi_rotation.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

double off_diagonl(double** array, int size){
    double sum = 0;
    for (int i=0;i<size;i++){
        for (int j=0;j<size;j++){
            if (j==i){
                continue;
            }
            if (array[i][j] != array[j][i]){
                throw std::runtime_error("Array not symmetric!");
            }
            sum = sum + array[i][j]*array[i][j];
        }
    }
    return std::sqrt(sum);
};

void rotate(double** array, double** eigvectors, int size, int l, int k, double c, double s){
    double akk,all,akl,aik,ail,eig_ik,eig_il;

    akk = array[k][k];
    all = array[l][l];
    akl = array[k][l];

    array[k][k] = akk*c*c  - 2*akl*c*s + all*s*s;
    array[l][l] = all*c*c +  2*akl*c*s + akk*s*s;
    array[k][l] = 0;
    array[l][k] = 0;

    for (int i=0;i<size;i++){
        eig_ik = eigvectors[i][k];
        eig_il = eigvectors[i][l];

        eigvectors[i][k] = c*eig_ik - s*eig_il;
        eigvectors[i][l] = c*eig_il + s*eig_ik;

        if (i==k or i==l){
            continue;
        }
        aik = array[i][k];
        ail = array[i][l];

        array[i][k] = aik*c - ail*s;
        array[k][i] = array[i][k];
        array[i][l] = ail*c + aik*s;
        array[l][i] = array[i][l];
    }
};




int jacobi_rotate(double** array, double* eigvals, double** eigvectors, int size, int maxiter, double tolerance){
    /*

    Inputs:
    double** array: a 2d array of values of size nxn
    double* eigavls: an array of length n to store eigenvalues in
    double** eigvectors: an identity matrix of the same shape as array. Will transform to contain the eigenvectors in its columns.
    int size: size n of the array
    double tolerance: the tolerance factor in the jacobi jacobi_rotation

    Output:
    int: the number of iterations that was needed
    */
    int number_of_iterations = 0;
    int k = 1;
    int l = 0;
    double tau,t,c,s;
    //double tau = 0;
    //double t = 0;
    //double c = 0;
    //double s = 0;

    while (off_diagonl(array,size)>=tolerance && maxiter>number_of_iterations){
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
            t = 1.0/(tau + std::sqrt(1 + tau*tau));
        }else{
            t = -1.0/(-tau + std::sqrt(1 + tau*tau));
        }
        c = 1/std::sqrt(1+t*t);
        s = t*c;
        }else{
            c = 1;
            s = 0;
        }

        //rotate
        rotate(array,eigvectors,size,l,k,c,s);
        number_of_iterations++;

        /* uncomment to print the rotations performed
        std::cout << "rotation\n";
        std::cout << std::fixed;
        std::cout << std::setprecision(2);
        for (int i = 0;i<size;i++){
            for (int j = 0;j<size;j++){
                std::cout << array[i][j] << " ";
            }
            std::cout << "\n";
        }
        */
    }
    for (int i=0;i<size;i++){
    eigvals[i] = array[i][i];
    }
return number_of_iterations;
};



void write_to_file(std::string filename, double** eigvectors, double* eigvals, int iters, int size){
    /*
    writes to filename
    */

    std::ofstream outfile;

    outfile.open(filename);

    outfile << iters << "\n";
    for (int i=0;i<size+1;i++){
        for (int j=0; j<size;j++){
            if (i==0){
                outfile << eigvals[j] << " ";
                continue;
            }
            outfile << eigvectors[j][i-1] << " ";

        }
        outfile << "\n";
    }
    outfile.close();
};

