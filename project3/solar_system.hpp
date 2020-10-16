#include<string>
#include<armadillo>
#ifndef SOLAR_SYSTEM_H
#define SOLAR_SYSTEM_H




void read_inital_condition(std::string filename,int body_index, double* mass, arma::Col<double>* initial_position, arma::Col<double>* initial_velocity);
void read_initial_conditions(std::string filename, int number_of_bodies, int body_indices[], arma::Col<double>* masses, arma::Mat<double>* initial_position, arma::Mat<double>* initial_velocity);

void Euler(int number_of_timesteps, double dt, int number_of_bodies, arma::Col<double> masses, arma::Mat<double>* initial_positions, arma::Mat<double>* initial_velocities, arma::Cube<double>* positions, arma::Cube<double>* velocities);
void VelocityVerlet(int number_of_timesteps, double dt, int number_of_bodies, arma::Col<double> masses, arma::Mat<double>* initial_positions, arma::Mat<double>* initial_velocities, arma::Cube<double>* positions, arma::Cube<double>* velocities);

arma::Mat<double> force(int number_of_bodies,arma::Col<double> masses, arma::Mat<double> positions);

void write_to_file(int dims, int number_of_bodies, int number_of_timesteps, arma::Cube<double> positions, std::string filename);

#endif