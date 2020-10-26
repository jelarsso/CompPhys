#include "solar_system.hpp"
#include<string>
#include<armadillo>

const double pi = 3.14159265358979;


int main(int argv, char *argc[]){
    /*
    This program simulates the orbit of Mercury with the General relativity correcational term.
    The command line arguments are, simulation_length (in years), dt (in years), x0, y0, vx0, vy0 (the initial position and velocity of the planet)
    */
    double simulation_length = std::atof(argc[1]);
    double dt = std::atof(argc[2]);

    double x0 = std::atof(argc[3]);
    double y0 = std::atof(argc[4]);
    double vx0 = std::atof(argc[5]);
    double vy0 = std::atof(argc[6]);
    

    int timesteps = (int)(simulation_length/dt);
    std::cout << timesteps << "\n";
    int *indices;
    indices = new int[1];
    indices[0] = 1;

    arma::Col<double> mass(1);
    mass(0) = 0.0553; 
    
    arma::Mat<double> initial_pos(2,1);
    arma::Mat<double> initial_vel(2,1);

    initial_pos(0,0)=x0;
    initial_pos(1,0)=y0;

    initial_vel(0,0)=vx0;
    initial_vel(1,0)=vy0;
    
    

    SolarSystem sol(2,1,indices,2,true,true);
    sol.set_initial_conditions(initial_pos,initial_vel,mass);
    sol.get_init();
    sol.VelocityVerletMercury(timesteps,dt);
    sol.write_to_file("output.data");
    return 0;
}