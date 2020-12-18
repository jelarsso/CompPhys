#include "algorithms.hpp"
#include <armadillo>
#include <iostream>
#include <cmath>
#include<iomanip>
const double pi = 3.1415926535897932;

double radiohalflife(double t){
    // Not in use in this file.
    return t;
}

double fn(double t, double x){
    /*
    double t: the time at which to find the solution.
    double x: the position at which to find the solution.

    The analytical solution of the diffusion equation with boundaries set to 0, and the initial condition given as -x.

    return, double: the solution at t,x 
    */

    double tol = 1e-14; // double precision
    double val = 0;
    int n = 1;
    double this_term;
    
    while (n<1000){
        this_term = 2*(std::sin(pi*n) - pi*n*std::cos(pi*n))/(pi*pi*n*n)*std::exp(-n*n*pi*pi*t)*std::sin(n*pi*x);
        val-= this_term;
        n++;
    }
    
    return val;
}


double init_func(double x, double y){
    //Not used in this file.
    return 0;
}

double init_func(double x){
    /*
    double x: the position at which to return the initial condition.
    
    return, double: the initial condition for the algorithms at x
    */
    return -x;
}

int main(int argc, char* argv[]){
    double dx = std::atof(argv[1]);
    double L = 1;
    int nt = std::atoi(argv[2]);

    int nx = (int) (L/dx);
    double alpha = std::atof(argv[3]);
    std::cout << dx << " " << nx << "\n";

    arma::Col<double> usolve(nx+1,arma::fill::zeros);
    forward_euler(nx,dx,alpha,nt,&usolve,"forward_euler.data");
    usolve.zeros();
    backward_euler(nx,dx,alpha,nt,&usolve,"backward_euler.data");
    usolve.zeros();
    cranky_nicholson(nx,dx,alpha,nt,&usolve,"cnicholson.data");
    close();

    std::cout << "nx = "<< nx+1 << " nt = " << nt+1 << std::endl;

    arma::Col<double> u(nx+1,arma::fill::zeros);
    for (double t=0;t<nt+1; t++){
        for (int i=0;i<nx+1;i++) u(i) = fn(t*alpha*dx*dx,i*dx);
        output(nx,&u,"analytical.data");
    }
    close();

    return 0;
}