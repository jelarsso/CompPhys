#include "solar_system.hpp"
#include<armadillo>
#include<fstream>
#include<string>


SolarSystem::SolarSystem(int nbodies, int dimensions){
    number_of_bodies = nbodies;
    dims = dimensions;


};


void SolarSystem::read_state(std::string filename){
    std::fstream input_file;


    initial_positions = new arma::Mat<double>(dims,number_of_bodies);
    initial_velocities = new arma::Mat<double>(dims,number_of_bodies);

    //loop through and read file
};


void SolarSystem::solve(int number_of_timesteps, int dt, bool method){

    positions = new arma::Cube<double>(number_of_timesteps,dims,number_of_bodies);
    velocities = new arma::Cube<double>(number_of_timesteps,dims,number_of_bodies);

    for (int i=0;i<number_of_timesteps;i++){
        positions(i) = positions(i-1) + force(positions);
    }    

};
