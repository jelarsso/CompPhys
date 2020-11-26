#include<armadillo>
#include<iostream>


double init_func(double);
void backward_euler(int, double, double, int);
void lusolve(double, double, double, int, arma::Col<double>*, arma::Col<double>*);




void lusolve(double a, double b, double c, int n, arma::Col<double>* rhs, arma::Col<double>* sol){
    // solution is in sol
    arma::Col<double> l(n),d(n+1), u(n+1);
    d(0) = b;
    for (int i = 1; i<n+1; i++){
        l(i-1) = a/d(i-1);
        d(i) = b - l(i-1)*c;
    }

    u(0) = (*rhs)(0);
    for (int i = 1; i<n+1; i++){
        u(i) = (*rhs)(i) - u(i-1)*l(i-1);
    }
    
    u(n) = u(n)/d(n);
    (*sol)(n)  = u(n);
    for (int i = n-1; i>=1; i--){
        (*sol)(i) = (u(i) - (*sol)(i+1)*c)/d(i);
    }
};

void forward_euler(int n, double dx, double alpha, int number_of_steps_t){
    arma::Col<double> u(n+1), usolve(n+1, arma::fill::zeros);
    u(0) = 0;
    u(n) = 0;
    for (int i=1; i<n; i++){
        u(i) = init_func(dx*i);
    }
    
    for (int t=0; t<number_of_steps_t;t++){
        for (int i=0; i>n; i++){
            usolve(i+1) = alpha*u(i-1) + (1-2*alpha)*u(i) + alpha*u(i+1);
        }
    }

};

void backward_euler(int n, double dx,double alpha, int number_of_steps_t){
    arma::Col<double> u(n+1),uprev(n+1);
    u(0) = 0;
    uprev(0) = 0;
    u(n) = 0;
    uprev(n) = 0;
    for (int i=1; i<n; i++){
        u(i) = init_func(dx*i);
        uprev(i) = init_func(dx*i);
    }
    
    for (int i = 1; i < number_of_steps_t; i++){
        lusolve(-alpha,1+2*alpha,-alpha, n+1, &uprev, &u);
        u(0) = 0;
        u(n) = 0;
        for (int j = 0; j<n+1;j++) uprev(j) = u(j);
    }
};


