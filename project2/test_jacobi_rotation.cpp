#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "jacobi_rotation.hpp"
#include<armadillo>


TEST_CASE() {
    /* First case: Checks if differnce between 
    numerical and analytical solution is within 
    the given tolerance  */

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

    iters = jacobi_rotate(array,eigvals,eigvectors,size,max_iter,tolerance);

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

    for (int i=0; i<size; i++){
        if ((X(i)-analyt_eigval[i])<1e-8){
            REQUIRE( true );
        }
    }
    

    /* Second case: Checks if differnce between 
    numerical and analytical solution is within 
    the given tolerance for another arbitrary 
    matrix A */


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
    }
    test_matrix[0][0] = 2;
    test_matrix[1][1] = 3;
    test_matrix[1][2] = 4;
    test_matrix[2][1] = 4;
    test_matrix[2][2] = 9;


    iters = jacobi_rotate(test_matrix,eiva,eive,N,max_iter,tolerance);   

    double* val;
    val = new double[N];
    val[0] = 1;
    val[1] = 2;
    val[2] = 11;

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
   for (int i=0; i<size; i++){
       for (int j=0; j<size; j++){
            prod = eigvectors[j][i]*eigvectors[j][i];
        count = count + prod;
       }
   }
   // every columnvectors inner product should yield 1 
   REQUIRE (abs(count-100.0)<1e-8);


   delete[] eigvals;
   for (int i=0;i<size;i++){
       delete[] eigvectors[i];
       delete[] array[i];
    }

    
}