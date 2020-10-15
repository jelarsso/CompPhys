#include "solar_system.hpp"
#include<string>
#include<armadillo>



int main(int argv, char *argc[]){
    int timesteps = std::atoi(argc[1]);
    double dt = std::atof(argc[2]);

    std::string filename = "initial_conditions.data";
    int number_of_bodies = 2;
    int dims = 3;
    int index_of_bodies[] = {0,1};
    
    arma::Mat<double> inital_pos(dims,number_of_bodies,arma::fill::zeros);
    arma::Mat<double> inital_vel(dims,number_of_bodies,arma::fill::zeros);

    arma::Cube<double> pos(dims,number_of_bodies,timesteps, arma::fill::zeros);
    arma::Cube<double> vel(dims,number_of_bodies,timesteps, arma::fill::zeros);
    arma::Col<double> masses = {330000,1};

    inital_pos(0,1) = 1;
    inital_vel(1,1) = 2*3.14159265358979*std::sqrt(masses(0)/1);

    //read_initial_conditions(filename,number_of_bodies,index_of_bodies,&inital_pos,&inital_vel);
    //VelocityVerlet(10, 0.01, 3, &masses, &inital_pos, &inital_vel, &pos, &vel);

    //force(number_of_bodies,masses,inital_pos);
    Euler(timesteps, dt, number_of_bodies, masses, &inital_pos, &inital_vel, &pos, &vel);
    write_to_file(dims,number_of_bodies,timesteps,pos,"output.data");

    return 0;
}