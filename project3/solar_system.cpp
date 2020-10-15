#include "solar_system.hpp"
#include<armadillo>
#include<fstream>
#include<string>
#include<cmath>

const double pi = 3.14159265358979;

void read_inital_condition(std::string filename,int body_index, arma::Col<double>* initial_position, arma::Col<double>* initial_velocity){
    std::ifstream input_file;
    input_file.open(filename);
    std::string line;
    if (input_file.is_open()){
        for(int i=0;i<=3*body_index;i++){
            getline(input_file,line);
        }
        getline(input_file,line);
        for (int i=0;i<3;i++){
            initial_position->at(i) = stod(line.substr(3+(i*26),3+(i*26)+21));
        }
        getline(input_file,line);
        for (int i=0;i<3;i++){
            initial_velocity->at(i) = stod(line.substr(3+(i*26),3+(i*26)+21))*365.249;
        }
    }
    input_file.close();
};



void read_initial_conditions(std::string filename, int number_of_bodies, int body_indices[], arma::Mat<double>* initial_positions, arma::Mat<double>* initial_velocities){
    arma::Col<double> initial_velocity(3);
    arma::Col<double> initial_position(3);
    for (int i=0;i<number_of_bodies;i++){
        read_inital_condition(filename, body_indices[i], &initial_position, &initial_velocity);
        initial_positions->col(i) = initial_position;
        initial_velocities->col(i) = initial_velocity;
    }
};

arma::Mat<double> force(int number_of_bodies,arma::Col<double> masses, arma::Mat<double> positions){
    int dims = 3;
    arma::Mat<double> forces(dims,number_of_bodies,arma::fill::zeros);
    arma::Cube<double> all_forces(dims,number_of_bodies,number_of_bodies,arma::fill::zeros);
    double G = 4*pi*pi;

    for (int objA=0;objA<number_of_bodies;objA++){
        for (int objB=0; objB<objA;objB++){
            double distAB = arma::norm(positions.col(objA) - positions.col(objB));
            all_forces.slice(objB).col(objA) = -G*masses(objA)*masses(objB)*(positions.col(objA) - positions.col(objB))/std::pow(distAB,3);
            all_forces.slice(objA).col(objB) = -all_forces.slice(objB).col(objA);
        }
    }
    forces = arma::sum(all_forces,2);
    for(int i=0;i<number_of_bodies;i++){
        forces.col(i)/=masses(i); //conv to acc
    }
    return forces;
};

void VelocityVerlet(int number_of_timesteps, double dt, int number_of_bodies, arma::Col<double> masses, arma::Mat<double>* initial_positions, arma::Mat<double>* initial_velocities, arma::Cube<double>* positions, arma::Cube<double>* velocities){
    positions->slice(0) = *initial_positions;
    velocities->slice(0) = *initial_velocities;
    arma::Mat<double> prev_force = force(number_of_bodies,masses,positions->slice(0));
    arma::Mat<double> this_force;

    for (int i=1;i<number_of_timesteps;i++){
        positions->slice(i) = positions->slice(i-1) + (velocities->slice(i-1))*dt + (0.5*dt*dt)*prev_force;
        positions->slice(i).col(0).zeros();//zero out sun
        this_force = force(number_of_bodies,masses,positions->slice(i));
        velocities->slice(i) = velocities->slice(i-1) + (0.5*dt)*(prev_force + this_force);
        velocities->slice(i).col(0).zeros();//zero out sun

        prev_force = this_force;
    }
};


void Euler(int number_of_timesteps, double dt, int number_of_bodies, arma::Col<double> masses, arma::Mat<double>* initial_positions, arma::Mat<double>* initial_velocities, arma::Cube<double>* positions, arma::Cube<double>* velocities){
    positions->slice(0) = *initial_positions;
    velocities->slice(0) = *initial_velocities;
    for (int i=1;i<number_of_timesteps;i++){
        positions->slice(i) = positions->slice(i-1) + dt*velocities->slice(i-1);
        velocities->slice(i) = velocities->slice(i-1) + dt*force(number_of_bodies,masses,positions->slice(i));
        velocities->slice(i).col(0).zeros();
    }
};

void write_to_file(int dims, int number_of_bodies, int number_of_timesteps, arma::Cube<double> positions, std::string filename){
    std::ofstream output_file;
    output_file.open(filename);

    for (int k = 0; k<number_of_timesteps; k++){
        for (int j = 0;j<number_of_bodies; j++){
            for (int i = 0;i<dims; i++){
                output_file << positions(i,j,k);
                output_file << " ";
            }
            output_file << " ";
        }
        output_file << "\n";
    }
    output_file.close();
};
