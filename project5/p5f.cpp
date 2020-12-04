#include "algorithms.hpp"
#include <armadillo>
#include <iostream>
#include <cmath>


double init_func(double x, double y){
    if ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < 0.1){
        return 1;
    }else{
        return 0;
    }
}

double init_func(double x){
    return -x;
}

int main(int argc, char* argv[]){
    double dx = std::atof(argv[1]);
    double L = 1;
    int nt = std::atoi(argv[2]);

    int nx = (int) (L/dx);
    double alpha = 0.25;
    std::cout << dx << " " << nx << "\n";

    arma::Mat<double> usolve(nx+1,nx+1,arma::fill::zeros);
    forward_euler2d(nx,dx,alpha,nt,&usolve,"fe2d.data");

    return 0;
}