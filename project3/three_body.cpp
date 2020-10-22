#include "solar_system.hpp"
#include<string>
#include<armadillo>



int main(int argv, char *argc[]){
    double simulation_length = std::atof(argc[1]);
    double dt = std::atof(argc[2]);
    int timesteps = (int)(simulation_length/dt)+1;
    int *indices;
    indices = new int[1];
    indices[0] = 3;
    indices[1] = 5;

    SolarSystem sol(3,2,indices,2,false,true);

    sol.VelocityVerlet(timesteps,dt);
    sol.write_to_file("output.data");
    return 0;
}