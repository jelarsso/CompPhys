#include "algorithms.hpp"
#include <armadillo>
#include <iostream>
#include <cmath>


double init_func(double x, double y){
    return 0;
}

double init_func(double x){
    return -x;
}

int main(int argc, char* argv[]){
    double dx = std::atof(argv[1]);
    double L = 1;
    int nt = std::atoi(argv[2]);

    int nx = (int) (L/dx);
    double alpha = 0.5;
    std::cout << dx << " " << nx << "\n";

    arma::Col<double> usolve(nx+1,arma::fill::zeros);
    forward_euler(nx,dx,alpha,nt,&usolve,"forward_euler.data");
    usolve.zeros();
    backward_euler(nx,dx,alpha,nt,&usolve,"backward_euler.data");
    usolve.zeros();
    cranky_nicholson(nx,dx,alpha,nt,&usolve,"cnicholson.data");

    return 0;
}