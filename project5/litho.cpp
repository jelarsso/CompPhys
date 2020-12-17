#include "algorithms.hpp"
#include <armadillo>
#include <iostream>
#include <cmath>



double length_scale = 150085.68980990254; //m
double time_scale = 1e9*365*24*60*60; //s
double temp_scale = 1; //C

double radiohalflife(double t){
    double q = (0.5e-6/2.5*length_scale*length_scale/temp_scale)*(0.4*std::pow(0.5,t/4.47) + 0.4*std::pow(0.5,t/14.0) + 0.2*std::pow(0.5,t/1.25));
    return q;
}


double init_func(double x, double y){
        return 0;
}

double init_func(double x){
    return -x;
}

int main(int argc, char* argv[]){
    //Physical variables
    double L = 120e3; //m
    double upper_layer_boundary = 20e3; //m
    double lower_layer_boundary = 40e3; //m
    double surface_boundary = 8; // C
    double core_boundary = 1300; // C
    double qprod_lower = 0.05e-6; // J/sm3
    double qprod_middle = 0.35e-6; // J/sm3
    double qprod_upper = 1.4e-6; // J/sm3
    double width_subduction = 150e3; // km

    //Simulation variables
    double dx = std::atof(argv[1]); // in 150km units
    int nx = (int) (L/length_scale/dx);
    int ny = nx*2;
    double sim_time = std::atof(argv[2]); // in Gy
    double alpha = 0.25;
    double dt = alpha*dx*dx;
    int nt = (int) (sim_time/dt);

    //Scaling of variables:
    double qupper,qlower,qmiddle;
    qlower = qprod_lower/2.5*length_scale*length_scale/temp_scale;
    qmiddle = qprod_middle/2.5*length_scale*length_scale/temp_scale;
    qupper = qprod_upper/2.5*length_scale*length_scale/temp_scale;
    double s_bound, c_bound;
    s_bound = surface_boundary/temp_scale;
    c_bound = core_boundary/temp_scale;
    double ulb,llb;
    ulb = upper_layer_boundary/length_scale;
    llb = lower_layer_boundary/length_scale;
    double s_width;
    s_width = width_subduction/length_scale;

    std::cout << "Simulating a grid " << nx+1 << " by " << ny+1 << " for T = " << sim_time << " Gy \n";
    std::cout << "nt = " << nt + 1 << "\n";
    std::cout << "dt = " << dt << " dx = " << dx << "\n\n";
    std::cout << "Scaled units: rhocp/k = " << 1 << " qk_upper = " << qupper << "\n";

    //setup for simulations
    arma::Mat<double> usolve(nx+1,ny+1,arma::fill::zeros);

    for (int i = 0;i<ny+1;i++){
        usolve(0,i) = s_bound;
        usolve(nx,i) = c_bound;
    }

    //simulate:
    forward_euler2d_litho_pb(nx, ny, dx, alpha,  ulb, llb, 0, 0, 0, nt, &usolve, "litho_no_Q.data");
    forward_euler2d_litho_pb(nx, ny, dx, alpha,  ulb, llb, qupper, qmiddle, qlower, nt, &usolve, "litho_Q_pb.data");
    forward_euler2d_litho(nx, ny, dx, alpha,  ulb, llb, qupper, qmiddle, qlower, s_width, nt, &usolve, "litho_enriched.data");


    return 0;
}



