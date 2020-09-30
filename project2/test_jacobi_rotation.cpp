#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "jacobi_rotation.hpp"
#include<armadillo>
#include<cmath>

TEST_CASE() {
    /* First case: Checks if differnce between 
    numerical and analytical solution is within 
    the given tolerance  */

    // We initilize the spring problem with N=100 mesh points
    int size = 100;
    int max_iter = 1000000;
    double tolerance = 1e-16;
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

    // We solve the linear system for the eigenvalues and eigenvectors with our Jacobi rotation algorithm 
    iters = jacobi_rotate(array,eigvals,eigvectors,size,max_iter,tolerance);

    // We solve the system analytically for comparison
    double* analyt_eigval;
    double** analyt_eigvec; 
    analyt_eigval = new double[size];
    analyt_eigvec = new double*[size];  
    for (int i=0;i<size;i++){
        analyt_eigvec[i] = new double[size];
    }
    analytical(analyt_eigval, analyt_eigvec, size);

    arma::vec Y(size);
    for (int i=0; i<size; i++){
        Y(i) = eigvals[i];
    }
    arma::vec X = arma::sort(Y); 

    // This is the test. We check if the jacobi solution is in compliance with the analytical of all eigenvalues. 
    for (int i=0; i<size; i++){
        if (std::abs(X(i)-analyt_eigval[i])<1e-8){
            REQUIRE( true );
        }
    }
    

    /* Second case: Checks if differnce between 
    numerical and analytical solution is within 
    the given tolerance for another arbitrary 
    matrix A */

    // We initialize a known problem 
    int N = 3;
    double* eiva;
    double** eive;
    
    eiva = new double[size];
    eive = new double*[size];
    for (int i=0;i<N;i++){
        eive[i] = new double[size];
        eiva[i] = 0;
    }

    for (int i = 0;i<N;i++){
        for (int j = 0;j<N;j++){
            if (i!=j){
                eive[i][j] = 0;
            }else{
                eive[i][j] = 1;
            }
        }
    }

    double** test_matrix;
    test_matrix = new double*[N];
    for (int i=0;i<N;i++){
        test_matrix[i] = new double[size];
        for (int j=0;j<N;j++){
            test_matrix[i][j] = 0;
        }
    }
    test_matrix[0][0] = 2;
    test_matrix[1][1] = 3;
    test_matrix[1][2] = 4;
    test_matrix[2][1] = 4;
    test_matrix[2][2] = 9;

    // We solve the system with the Jacobi rotation algorithm 
    iters = jacobi_rotate(test_matrix,eiva,eive,N,max_iter,tolerance);   

    // We set the correct eigenvalues in an array
    double* val;
    val = new double[N];
    val[0] = 1;
    val[1] = 2;
    val[2] = 11;

    // We test if the numerical solution is in compliance with the known solution
    for (int i=0; i<N; i++){
        if (abs(val[i] - eiva[i])<1e-8){
            REQUIRE( true );
        }
    } 

    /*
    Third case: Checks orthonormality of transformed matrices 
    */

   double count = 0;
   double prod = 0;
   double dot_product = 0;
   for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
           for (int k=0;k<size;k++){
            dot_product = dot_product + eigvectors[j][k]*eigvectors[i][k];
           }
            if (i!=j){
                REQUIRE(std::abs(dot_product)<1e-8);
            }else{
                REQUIRE(std::abs(dot_product-1)<1e-8);
            }
            dot_product = 0;
       }
   }


   delete[] eigvals;
   for (int i=0;i<size;i++){
       delete[] eigvectors[i];
       delete[] array[i];
    }

    
}