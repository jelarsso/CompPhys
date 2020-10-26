#include "solar_system.hpp"
#include<string>
#include<armadillo>



int main(int argv, char *argc[]){
    /*
    Simulates the Earth's orbit around the sun starting at (1,0)AU with velcoity (0,start_vel)AU/year.
    The command line arguments are, simulation_length (in years), dt (timestep-length in years), start_vel (the initial velcoity in the y direction), 
    beta (the parameter for modified gravity).
    */
    double simulation_length = std::atof(argc[1]);
    double dt = std::atof(argc[2]);
    double start_vel = std::atof(argc[3]);
    double beta = 2; // std::atof(argc[4]);
    int timesteps = (int)(simulation_length/dt)+1;
    int *indices;
    indices = new int[1];
    indices[0] = 1;

    arma::Mat<double> init_pos = arma::Mat<double>(2,1);
    arma::Mat<double> init_vel = arma::Mat<double>(2,1);
    arma::Col<double> masses = arma::Col<double>(1);
    init_pos(0,0) = 1;
    init_pos(1,0) = 0;
    init_vel(0,0) = 0;
    init_vel(1,0) = start_vel;
    masses(0) = 1;
    

    SolarSystem sol(2,1,indices,beta,false,true);
    sol.set_initial_conditions(init_pos,init_vel,masses);
    sol.VelocityVerlet(timesteps,dt);
    sol.write_to_file("output.data");

    return 0;
}