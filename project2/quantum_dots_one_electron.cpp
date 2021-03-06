#include "jacobi_rotation.hpp"
#include<iostream>
#include<string>
#include<iomanip>


int main(int argc, char *argv[]){
    /*
    This program solves the one electron problem
    */
    int size = std::atoi(argv[1]);
    double tol = std::atof(argv[2]);
    //double rho_max = std::atof(argv[3]); for testing different rho_maxes.
    int iters;
    

    double rho_max = 4.5; 

    double h = rho_max/(1.0*size);
    std::cout << "h = " << h << std::endl;
    double d = 2.0/(h*h);
    double a = -1.0/(h*h);

    double* eigvals;
    double** array;
    double** eigvectors;
    
    array = new double*[size];
    eigvals = new double[size];
    eigvectors = new double*[size];
    for (int i=0;i<size;i++){
        array[i] = new double[size];
        eigvectors[i] = new double[size];
        for (int j=0;j<size;j++){
            array[i][j]=0;
        }
    }

    for (int i = 0;i<size;i++){
        array[i][i] = d + (i+1)*(i+1)*h*h;
        if(i<size-1){
        array[i][i+1] = a;
        array[i+1][i] = a;
        }
        for (int j = 0;j<size;j++){
            if (i!=j){
                eigvectors[i][j] = 0;
            }else{
                eigvectors[i][j] = 1;
            }            
        }
    }

    iters = jacobi_rotate(array,eigvals,eigvectors,size,10000000,tol);

    std::cout << iters << std::endl;

    std::string filename = "eigen.data";
    write_to_file(filename,eigvectors,eigvals,iters,size);

    for (int i=0;i<size;i++){
        std::cout << std::setprecision(15) << eigvals[i] << std::endl;
    }

    delete[] eigvals;
    for (int i=0;i<size;i++){
        delete[] eigvectors[i];
        delete[] array[i];
    }
        
    return 0;
};