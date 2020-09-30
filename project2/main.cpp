#include "jacobi_rotation.hpp"
#include<iostream>
#include<string>



int main(int argc, char *argv[]){
    /*
    This program calculates the simple spring with fixed ends, prints the number of iterations, and writes the solution to file.
    */
    int size = std::atoi(argv[1]);
    double tol = std::atof(argv[2]);
    int iters;
    
    double h = 1/(1.0*size);
    double d = 2/(h*h);
    double a = -1/(h*h);

    double* eigvals;
    double** array;
    double** eigvectors;
    
    array = new double*[size];
    eigvals = new double[size];
    eigvectors = new double*[size];
    for (int i=0;i<size;i++){
        array[i] = new double[size];
        eigvectors[i] = new double[size];
    }

    for (int i = 0;i<size;i++){
        for (int j = 0;j<size;j++){
            array[i][j] = 0;
            if (i!=j){
                eigvectors[i][j] = 0;
            }else{
                eigvectors[i][j] = 1;
            }            
        }
        array[i][i] = d;
        if(i>0){
        array[i][i-1] = a;
        array[i-1][i] = a;
        }
    }

    iters = jacobi_rotate(array,eigvals,eigvectors,size,10000000,tol);

    std::cout << iters << std::endl;

    for (int i=0; i<size; i++){
        std::cout << eigvals[i] << "\n";
    }

    std::string filename = "data.out";
    write_to_file(filename,eigvectors,eigvals,iters,size);

    delete[] eigvals;
    for (int i=0;i<size;i++){
        delete[] eigvectors[i];
        delete[] array[i];
    }
        
    return 0;
};