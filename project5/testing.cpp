#include "algorithms.hpp"
#include <armadillo>
#include <iostream>
#include <cmath>


double init_func(double x, double y){
    return 0;
}

double init_func(double x){
    return 0;
}

int main(){
    int nx = 100;
    double dx = 0.1;
    double alpha = 0.001;
    int nt = 100000;
    arma::Col<double> usolve(nx+1,arma::fill::zeros);

    forward_euler(nx,dx,alpha,0,1,nt,&usolve,"dump.data");

    return 0;
}