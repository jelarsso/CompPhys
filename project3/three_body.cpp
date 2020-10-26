#include "solar_system.hpp"
#include<string>
#include<armadillo>



int main(int argv, char *argc[]){
    double simulation_length = std::atof(argc[1]);
    double dt = std::atof(argc[2]);
    int timesteps = (int)(simulation_length/dt)+1;
    int *indices;
    indices = new int[2];
    indices[0] = 3;
    indices[1] = 5;
    indices[2] = 0;

    SolarSystem sol(3,3,indices,2,false,false);

    sol.change_reference_to_cm();
    sol.VelocityVerlet(timesteps,dt);
    sol.write_to_file("output.data");
    return 0;
}