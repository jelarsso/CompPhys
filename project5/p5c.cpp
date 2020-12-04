#include "algorithms.hpp"
#include <armadillo>
#include <iostream>
#include <cmath>

const double pi = 3.1415926535897932;

double radiohalflife(double t){
    return t;
}

double fn(double t, double x){
    int terms = 1000;
    double val = 0;
    for (int n=1;n<terms;n++){
        val-= 2*(std::sin(pi*n) - pi*n*std::cos(pi*n))/(pi*pi*n*n)*std::exp(-n*n*pi*pi*t)*std::sin(n*pi*x);
    }
    return val;
}


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
    close();

    std::cout << nx+1 << " " << nt+1 << std::endl;

    arma::Col<double> u(nx+1,arma::fill::zeros);
    for (double t=0;t<nt+1; t++){
        for (int i=0;i<nx+1;i++) u(i) = fn(t*alpha*dx*dx,i*dx);
        output(nx,&u,"analytical.data");
    }
    u.print();
    close();

    return 0;
}