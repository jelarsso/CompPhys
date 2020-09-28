#include "jacobi_rotation.hpp"
#include<iostream>
#include<string>



int main(int argc, char *argv[]){
    int size = std::atoi(argv[1]);
    double tol = std::atof(argv[2]);
    double omega_r = std::atof(argv[3]);
    int iters;
    

    double rho_max = 10; //size ser ut til å måtte være mye større enn rhomax for å gi løsningen

    double h = rho_max/(1.0*size);
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
        array[i][i] = d + omega_r*omega_r*(i+1)*(i+1)*h*h + 1/(1.0*(i+1)*h);
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

    iters = jacobi_rotate(array,eigvals,eigvectors,size,10000000,1e-16);


    std::string filename = "two_electron.data";
    write_to_file(filename,eigvectors,eigvals,iters,size);


    delete[] eigvals;
    for (int i=0;i<size;i++){
        delete[] eigvectors[i];
        delete[] array[i];
    }
        
    return 0;
};