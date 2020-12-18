#include "algorithms.hpp"
#include <armadillo>
#include <iostream>
#include <cmath>

const double pi = 3.1415926535897932;


double init_func(double x, double y){
    /*if ((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < 0.1){
        return 1;
    }else{
        return 0;
    }*/
    return -x;
}

double init_func(double x){
    return -x;
}

double radiohalflife(double t){
    return t;
}


double coeffs(int nx, int ny, double x, double y){
    double lx = 1;
    double ly = 1;
    return (4/(lx*ly))*((ly/(ny*pi))*std::cos(ny*pi)-ly/(ny*pi))*(-lx*lx/(pi*nx)*std::cos(nx*pi))*std::sin(nx*pi*x/lx)*std::sin(ny*pi*y/ly);
}

double fn(double t, double x, double y){
    double tol = 1e-14; // double precision
    double val = 0;

    double lx = 1;
    double ly = 1;
    
    for (int i=1;i<50;i++){
        for (int j=1; j<50; j++){
            val += coeffs(i,j,x,y)*std::exp(-i*i*pi*pi*t/lx/lx-j*j*pi*pi*t/ly/ly);
        }
    }
    return val;
}



int main(int argc, char* argv[]){
    double dx = std::atof(argv[1]);
    double L = 1;
    int nx = (int) (L/dx);
    int nt = std::atoi(argv[2]);
    double alpha = std::atof(argv[3]);

    std::cout << "dx = " << dx << "  nx = " << nx << " nt = " << nt << "\n";

    arma::Mat<double> usolve(nx+1,nx+1,arma::fill::zeros);
    forward_euler2d(nx,dx,alpha,nt,&usolve,"fe2d.data");

    usolve.zeros();
    for (int t=0; t<nt+1; t++){
        for (int ix=0; ix<nx+1;ix++){
            for (int iy=0; iy<nx+1;iy++){
                usolve(ix,iy) = fn(t*alpha*dx*dx,ix*dx,iy*dx);
            }
        }
        output(nx,nx,&usolve,"analytical2d.data");
    }
    close();

    return 0;
}