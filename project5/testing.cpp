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
    /*
    int nx = 100;
    double dx = 0.1;
    double alpha = 0.001;
    int nt = 100000;
    arma::Col<double> usolve(nx+1,arma::fill::zeros);

    forward_euler(nx,dx,alpha,0,1,nt,&usolve,"dump.data");
    */

    int nx = 100;
    double dx = .1;
    double alpha = .0001;
    int nt = 10000;
    arma::Col<double> usolve(nx+1,arma::fill::zeros);
    
    forward_euler(nx,dx,alpha,0,1,nt,&usolve,"fe2.data");
    backward_euler(nx,dx,alpha,0,1,nt,&usolve,"be2.data");
    cranky_nicholson(nx,dx,alpha,0,1,nt,&usolve,"cn2.data");

    return 0;
}