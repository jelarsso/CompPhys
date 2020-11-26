#include<armadillo>
#include<iostream>


double init_func(double);
void backward_euler(int, double, double, int);
void lusolve();



void forward_euler(){
    
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
    
    for (int i = 1; i < number_of_steps_t, i++){
        lusolve();
        u(0) = 0;
        u(n) = 0;
        for (int j = 0; j<n+1;j++) uprev(j) = u(j);
    }
};


