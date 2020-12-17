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


double fn(double t, double x, double y){
    double tol = 1e-14; // double precision
    double val = 0;
    int nx = 1;
    int ny = 1;
    double this_term;


    double lx = 1;
    double ly = 1;
    
    while (nx<100){
        while (ny<100){
            this_term = std::sin(nx*pi*x/lx)*std::sin(ny*pi*y/ly)*((-lx*lx/(nx*pi)*std::cos(nx*pi))*(-ly/(pi*ny)+ly/(ny*pi)*std::cos(pi*ny)));
            val -= 4/lx/ly*(this_term)*std::exp( - nx*nx*pi*pi*t/lx/lx  - ny*ny*pi*pi*t/ly/ly );
            ny++;
        }
        nx++;
    }
    
    return val;
}



int main(int argc, char* argv[]){
    double dx = std::atof(argv[1]);
    double L = 1;
    int nt = std::atoi(argv[2]);

    int nx = (int) (L/dx);
    double alpha = 0.25;
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