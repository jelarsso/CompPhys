#include "solar_system.hpp"
#include<armadillo>
#include<fstream>
#include<string>
#include<cmath>

const double pi = 3.14159265358979;
const double sun_mass = 332946.0487; // in units of earth masses
const double G1 = 4*pi*pi;
const double G = 4*pi*pi/sun_mass;

SolarSystem::SolarSystem(int dim, int nbodies, int *bindices, double beta_value, bool set_initial_conditions_manually, bool stationary_s){
    dims = dim;
    number_of_bodies = nbodies;
    body_indices = bindices;
    

    beta = beta_value;
    stationary_sun = stationary_s;

    initial_positions = arma::Mat<double>(dims,number_of_bodies, arma::fill::zeros);
    initial_velocities = arma::Mat<double>(dims,number_of_bodies,arma::fill::zeros);
    masses = arma::Col<double>(number_of_bodies,arma::fill::zeros);

    if (set_initial_conditions_manually == false){
        SolarSystem::read_initial_conditions("initial_conditions.data");
    }
};


void SolarSystem::set_initial_conditions(arma::Mat<double> start_positions, arma::Mat<double> start_velocities, arma::Col<double> mass){
    initial_positions = start_positions;
    initial_velocities = start_velocities;
    masses = mass;
};

void SolarSystem::read_initial_condition(std::string filename,int body_index, double* mass, arma::Col<double>* initial_position, arma::Col<double>* initial_velocity){
    std::ifstream input_file;
    input_file.open(filename);
    std::string line;
    if (input_file.is_open()){
        for(int i=0;i<=4*body_index;i++){
            getline(input_file,line);
        }
        getline(input_file,line);
        *mass = stod(line.substr(4,14));
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

void SolarSystem::read_initial_conditions(std::string filename){
    arma::Col<double> initial_velocity(dims);
    arma::Col<double> initial_position(dims);
    double mass;
    for (int i=0;i<number_of_bodies;i++){
        SolarSystem::read_initial_condition(filename, body_indices[i], &mass, &initial_position, &initial_velocity);
        initial_positions.col(i) = initial_position;
        initial_velocities.col(i) = initial_velocity;
        masses(i) = mass;
    }
};

arma::Mat<double> SolarSystem::force(arma::Mat<double> positions){
    arma::Mat<double> forces(dims,number_of_bodies,arma::fill::zeros);
    arma::Cube<double> all_forces(dims,number_of_bodies,number_of_bodies,arma::fill::zeros);
    
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

arma::Mat<double> SolarSystem::stat_sun_force(arma::Mat<double> positions){
    arma::Mat<double> forces(dims,number_of_bodies,arma::fill::zeros);
    for (int objA=0;objA<number_of_bodies;objA++){
        double distA = arma::norm(positions.col(objA));
        forces.col(objA) = -G1*(positions.col(objA))/std::pow(distA,(1+beta));
    }
    if (number_of_bodies>1){
        forces = forces + force(positions);
    }
    return forces;
};

void SolarSystem::VelocityVerlet(int number_of_timesteps, double dt_length){
    dt = dt_length;
    timesteps = number_of_timesteps;
    positions = arma::Cube<double> (dims,number_of_bodies,timesteps,arma::fill::zeros);
    velocities = arma::Cube<double> (dims,number_of_bodies,timesteps,arma::fill::zeros);
    positions.slice(0) = initial_positions;
    velocities.slice(0) = initial_velocities;
    arma::Mat<double> prev_force;
    arma::Mat<double> this_force;

    if (stationary_sun==false){
    prev_force = force(positions.slice(0));
    for (int i=1;i<number_of_timesteps;i++){
        positions.slice(i) = positions.slice(i-1) + (velocities.slice(i-1))*dt + (0.5*dt*dt)*prev_force;
        this_force = force(positions.slice(i));
        velocities.slice(i) = velocities.slice(i-1) + (0.5*dt)*(prev_force + this_force);
        prev_force = this_force;
    }}else{
    prev_force = stat_sun_force(positions.slice(0));
    for (int i=1;i<number_of_timesteps;i++){
        positions.slice(i) = positions.slice(i-1) + (velocities.slice(i-1))*dt + (0.5*dt*dt)*prev_force;
        this_force = stat_sun_force(positions.slice(i));
        velocities.slice(i) = velocities.slice(i-1) + (0.5*dt)*(prev_force + this_force);
        prev_force = this_force;
    }}
};


void SolarSystem::Euler(int number_of_timesteps, double dt_length){
    dt = dt_length;
    timesteps = number_of_timesteps;
    positions = arma::Cube<double>(dims,number_of_bodies,timesteps,arma::fill::zeros);
    velocities = arma::Cube<double>(dims,number_of_bodies,timesteps,arma::fill::zeros);
    positions.slice(0) = initial_positions;
    velocities.slice(0) = initial_velocities;
    if (stationary_sun==false){
    for (int i=1;i<number_of_timesteps;i++){
        positions.slice(i) = positions.slice(i-1) + dt*velocities.slice(i-1);
        velocities.slice(i) = velocities.slice(i-1) + dt*force(positions.slice(i-1));
    }}else{
    for (int i=1;i<number_of_timesteps;i++){
        positions.slice(i) = positions.slice(i-1) + dt*velocities.slice(i-1);
        velocities.slice(i) = velocities.slice(i-1) + dt*stat_sun_force(positions.slice(i));
    }}
};

void SolarSystem::write_to_file(std::string filename){
    std::ofstream output_file;
    output_file.open(filename);

    output_file << "# timesteps: " << timesteps << " dt: " << dt << " nbodies: " << number_of_bodies << " dims: " << dims << "\n";
    output_file << "#each column corresponds to for x11 v11 x21 v21 (x31) (v31) x12 v12 x22 v22 (x32) v(32) ... for each body\n"; 
    for (int k = 0; k<timesteps; k++){
        for (int j = 0;j<number_of_bodies; j++){
            for (int i = 0;i<dims; i++){
                output_file << positions(i,j,k);
                output_file << " ";
                output_file << velocities(i,j,k);
                output_file << " ";
            }
            output_file << " ";
        }
        output_file << "\n";
    }
    output_file.close();
};

void SolarSystem::change_reference_to_cm(){
    arma::Col<double> vel_cm = initial_velocities*masses /arma::sum(masses);
    for (int i=0; i<number_of_bodies; i++){
        initial_velocities.col(i) -= vel_cm;
    };
    
};


void SolarSystem::get_init(){
    initial_positions.print();
    std::cout << " " << std::endl;
    initial_velocities.print();
    std::cout << "\n";
    masses.print();
};

SolarSystem::~SolarSystem(){

};


/*
// FUNCTION IMPLEMENTATION OF 3A
void read_inital_condition(std::string filename,int body_index, double* mass, arma::Col<double>* initial_position, arma::Col<double>* initial_velocity){
    std::ifstream input_file;
    input_file.open(filename);
    std::string line;
    if (input_file.is_open()){
        for(int i=0;i<=4*body_index;i++){
            getline(input_file,line);
        }
        getline(input_file,line);
        *mass = stod(line.substr(4,14));
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



void read_initial_conditions(std::string filename, int number_of_bodies, int body_indices[], arma::Col<double>* masses, arma::Mat<double>* initial_positions, arma::Mat<double>* initial_velocities){
    arma::Col<double> initial_velocity(3);
    arma::Col<double> initial_position(3);
    double mass;
    for (int i=0;i<number_of_bodies;i++){
        read_inital_condition(filename, body_indices[i], &mass, &initial_position, &initial_velocity);
        initial_positions->col(i) = initial_position;
        initial_velocities->col(i) = initial_velocity;
        masses->at(i) = mass;
    }
};

arma::Mat<double> force(int number_of_bodies,arma::Col<double> masses, arma::Mat<double> positions){
    int dims = 3;
    arma::Mat<double> forces(dims,number_of_bodies,arma::fill::zeros);
    arma::Cube<double> all_forces(dims,number_of_bodies,number_of_bodies,arma::fill::zeros);
    

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
        //positions->slice(i).col(0).zeros();//zero out sun
        this_force = force(number_of_bodies,masses,positions->slice(i));
        velocities->slice(i) = velocities->slice(i-1) + (0.5*dt)*(prev_force + this_force);
        //velocities->slice(i).col(0).zeros();//zero out sun

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
*/