#include "jacobi_rotation.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>

double off_diagonl(double** array, int size){
    /*
    double** array, the array to the find the squared sum of the offdiagonal elements
    int size, the size of the aforementioned array

        Implementation of the function off(A) defined in the article.

    outputs:
    double: the value of the square root of the sum of the squared off diagonal elements.
    */
    double sum = 0;
    for (int i=0;i<size;i++){
        for (int j=0;j<size;j++){
            if (j==i){
                continue;
            }
            // can be uncommented to make sure that the array stays symmetric, does however not seem to be compatible with catch so it is commented.
            //if (array[i][j] != array[j][i]){
            //    throw std::runtime_error("Array not symmetric!");
            //}
            sum = sum + array[i][j]*array[i][j];
        }
    }
    return std::sqrt(sum);
};

void rotate(double** array, double** eigvectors, int size, int l, int k, double c, double s){
    /*
    Inputs:
    double** array, the array to perform one single rotation on. a 2d array of size equal to the parameter size
    double** eigvectors, the 2d array representing the eigenvectors of the same shape as array
    int size, the size of the array
    int l,k: the array element to perform the rotation around
    double c,s: the amount to rotate (c = cos theta and s = sin theta, where theta is the rotation angle)

        Performs a single rotation of the matrix around element (l,k).
    
    outputs:
    (none)
    */
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

        Perform the Jacobi algorithm on the matrix stored in array. The eigenvalues are found in the vector eigvals and eigenvectors are in eigvectors.
        maxiter is the max number of iterations that will be performed.
        tolerance is the tolerance of the jacobi algorithm. i.e. the number that is compared with the function off(A) defined in the article.

    Output:
    int: the number of iterations that was needed
    */
    int number_of_iterations = 0;
    int k = 1;
    int l = 0;
    double tau,t,c,s;

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
    string filename: a string of the filename to store the data in
    the other parameters are described in jacobi_rotate
    int iters: the number of iterations from the jacobi rotation.

        writes results to filename. The first line is the number of iterations, the second the eigenvalues, then the next lines are the eigenvectors.
    
    outputs:
    (none)
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


void analytical(double* analyt_eigval, double** analyt_eigvec, int N) {
    /*
    Solves analytically the spring problem

    Inputs:
    double* analyt_eigval: an array of length n to store eigenvalues in 
    double** analyt_eigvec: an identity matrix of the same shape as array. Will transform to contain the eigenvectors in its columns
    int N: size n of the array
    
    Outputs:
    */
    double h = 1.0/N;
    double d = 2.0/(h*h);
    double a = -1.0/(h*h);

    for (int i=1; i<(N+1); i++){
        analyt_eigval[i-1] = d + 2*a*std::cos(i*M_PI/(N+1));
        for (int j=1; j<(N+1); j++){
            analyt_eigvec[i-1][j-1] = std::sin(i*j*M_PI/(N+1));
        }
    }
};
