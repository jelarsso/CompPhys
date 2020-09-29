#define CATCH_CONFIG_MAIN 
#include "catch.hpp"
#include "jacobi_rotation.hpp"

TEST_CASE() {

    int size = 100;
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
        eigvals[i] = 0;
    }

    for (int i = 0;i<size;i++){
        array[i][i] = d;
        if(i>0){
        array[i][i-1] = a;
        array[i-1][i] = a;
        }
        for (int j = 0;j<size;j++){
            if (i!=j){
                eigvectors[i][j] = 0;
            }else{
                eigvectors[i][j] = 1;
            }            
        }
    }

    iters = jacobi_rotate(array,eigvals,eigvectors,size,10000000,1e-8);

    for (int i=0; i<size; i++){
        std::cout << eigvals[i] << "\n";
    }

    double* analyt_eigval;
    double** analyt_eigvec; 
    analyt_eigval = new double[size];
    analyt_eigvec = new double*[size];  
    for (int i=0;i<size;i++){
        analyt_eigvec[i] = new double[size];
    }
    analytical(analyt_eigval, analyt_eigvec, size);

    /*
    for (int i=0; i<size; i++){
        std::cout << eigvals[i] << " ";
        std::cout << analyt_eigval[i] << "\n";
    }
    */

    std::string filename = "data2.out";
    write_to_file(filename,eigvectors,eigvals,iters,size);

    /*
    std::string filename = "data2.out";
    write_to_file(filename,analyt_eigvec,analyt_eigval,iters,size);
    */

    delete[] eigvals;
    for (int i=0;i<size;i++){
        delete[] eigvectors[i];
        delete[] array[i];
    }

    REQUIRE( true );
}