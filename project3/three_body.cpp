#include "solar_system.hpp"
#include<string>
#include<armadillo>

const double pi = 3.14159265358979;


int main(int argv, char *argc[]){
    /*
    Simulates the orbit of Earth and Jupiter with the interaction between them. Reads the inital conditions from the initial_conditions.data file.
    The command line arguments are simulation_length (in years), dt (in years), jupiter_mass_factor (the factor to multply jupiters mass with)
    */
    double simulation_length = std::atof(argc[1]);
    double dt = std::atof(argc[2]);
    double jupiter_mass_factor = std::atof(argc[3]);
    int timesteps = (int)(simulation_length/dt);
    int *indices;
    indices = new int[2];
    indices[0] = 3;
    indices[1] = 5;

    arma::Col<double> mass(2);
    mass(0) = 1;
    mass(1) = 317.8*jupiter_mass_factor;
        

    SolarSystem sol(2,2,indices,2,false,true);
    sol.set_initial_conditions(mass);
    sol.get_init();
    std::cout << timesteps << " " << dt << std::endl;
    sol.VelocityVerlet(timesteps,dt);
    sol.write_to_file("output.data");
    return 0;
}