#include "algorithms.hpp"
#include <armadillo>
#include <iostream>
#include <cmath>


double init_func(double x){
    if (x>0.45 && x<0.55){
        return 1;
    }
    return 0;
}

int main(){
    int nx = 100;
    double dx = 0.01;
    double alpha = 0.001;
    int nt = 100000;
    arma::Col<double> usolve(nx+1);

    //forward_euler(nx,dx,alpha,nt,&usolve,"dump.data");
    cranky_nicholson(nx,dx,alpha,nt,&usolve,"dump.data");


    return 0;
}