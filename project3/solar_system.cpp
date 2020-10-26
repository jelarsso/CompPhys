#include "solar_system.hpp"
#include<armadillo>
#include<fstream>
#include<string>
#include<iomanip>
#include<cmath>

const double pi = 3.14159265358979;
const double sun_mass = 332946.0487; // in units of earth masses
const double G1 = 4*pi*pi;
const double G = 4*pi*pi/sun_mass;
const double c = 63239.7263;

SolarSystem::SolarSystem(int dim, int nbodies, int *bindices, double beta_value, bool set_initial_conditions_manually, bool stationary_s){
    /*
    Inputs:
    int dim: the number of dimensions to set up the simulations for. (2 or 3)
    int nbodies: the number of bodies to simulate the system with, a positive integer.
    int *bindices: pointer to the head of an array of size nbodies. The entries correspond to the indices in the initial_conditions.data file:
                    0: Sun, 1: Mercury ..... 8: Uranus, 9: Pluto
    double beta_value: set to 2 to use the normal inverse square law of gravity (only supported when stationary_s==true and using the Verlet method)
    bool set_initial_conditions_manually: toogles whether to read the initial conditions from the file initial_conditions.data. If true, the method SolarSystem::set_initial_conditions() must be used.
    bool stationary_s: If stationary_s == true, the sun's position is not integrated and assumed to be at (0,0,0) for every timestep. The index correpsonding to the Sun (0 in the initial_conditions.data) must not be used.

        Constructor for the class SolarSystem.

        Class SolarSystem:
            Description:
                The SolarSystem class is a class that handles the simulation of N-body systems. So far only the graviational force has been implemented,
                but adding other forces is as simple as creating a new method for the class and adding it to the integrators.
                The class has the ability to integrate one-body and two-body systems where the sun is a stationary particle at (0,0)
                After the class has been initialized either by manually reading setting the initial conditions or by simply reading them from a file,
                the system can be solved using either the Euler or VelocityVerlet method.
                Finally the results can be easily written to a specified file.
                Example usages can be found in earth_sun_system.cpp or three_body.cpp
            
            Attributes:
                int dims, the number of dimensions of the system.
                int timesteps, the number of timesteps that was simulated.
                double simulation_length, how long (in years) that the system was simulated.
                double dt, the timestep-length that was used in years.
                double beta, the parameter for the modified gravity.
                int number_of_bodies, how many bodies that were simulated
                int* body_indices, the indices of the bodies corresponding to the planets in the data file initial_conditions.data
                bool stationary_sun, whether the sun is kept stationary at (0,0,0) or not

                arma::Mat<double> initial_positions, an armadillo Matrix of size dims * number_of_bodies containg the initial position vectors
                arma::Mat<double> initial_velocities, an armadillo Matrix of size dims * number_of_bodies containg the intitial velocity vectors
                arma::Col<double> masses, a vector of length number_of_bodies that contains the masses of the particles in units of earth masses.

                arma::Cube<double> positions, an armadillo Cube matrix of size dims * number_of_bodies * timesteps that contains the solved position vectors for each body at each timestep.
                arma::Cube<double> velocities, an armadillo Cube matrix of size dims * number_of_bodies * timesteps that contains the solved velocity vectors for each body at each timestep.

    */

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
    /*
    Inputs:
    arma::Mat<double> start_positions, an armadillo Matrix of size dims * number_of_bodies containg the initial position vectors
    arma::Mat<double> start_velocities, an armadillo Matrix of size dims * number_of_bodies containg the intitial velocity vectors
    arma::Col<double> mass, a vector of length number_of_bodies that contains the masses of the particles in units of earth masses.

        Sets the intitial conditions of the system manually from the arrays passed as arguments.
    */
    initial_positions = start_positions;
    initial_velocities = start_velocities;
    masses = mass;
};

void SolarSystem::set_initial_conditions(arma::Col<double> mass){
    /*
    Inputs:
    arma::Col<double> mass, a vector of length number_of_bodies that contains the masses of the particles in units of earth masses.
        
        Overloaded method, only sets the masses of the system. Allows to override the masses that was read from initial_conditions.data
    */
    masses = mass;
};


void SolarSystem::read_initial_condition(std::string filename,int body_index, double* mass, arma::Col<double>* initial_position, arma::Col<double>* initial_velocity){
    /*
    std::string filename, the file to read from
    int body_index, which object to read from the file
    double *mass, a pointer to where the mass is stored
    arma::Col<double>* initial_position, initial_velocity, pointers to the armadillo vectors of length dims, where the inital conditions for object is loaded into.

        A private method not really meant to be used explicitly. This reads from a file with the same format as initial_conditions.data a single objects initial conditions.
    */
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
    /*
    std::string filename, the name of the file to read the initial conditions of the system from.

        Reads in the initial conditions of the bodies that were specified in the constructor.
    */
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
    /*
    arma::Mat<double> positions, a matrix of size dims * number_of_bodies containg the position vectors as column vectors.

        Calculates the all the N-body interactions between the particles, using the positions specified as the parameter.
    
    returns:
    arma::Mat<double> positions, a matrix of size dims * number_of_bodies containg the acceleration vectors as column vectors.
    */
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
    /*
    arma::Mat<double> positions, a matrix of size dims * number_of_bodies containg the position vectors as column vectors.

        Calculates the acceleration from the sun on all the bodies, and all the N-body interactions between the particles (except the sun), using the positions specified as the parameter.
    
    returns:
    arma::Mat<double> positions, a matrix of size dims * number_of_bodies containg the acceleration vectors as column vectors.
    */
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

arma::Mat<double> SolarSystem::mercury_force(arma::Mat<double> positions, arma::Mat<double> velocities){
    /*
    arma::Mat<double> positions, a matrix of size dims * number_of_bodies containg the position vectors as column vectors.

        Calculates the acceleration on a single object including the general relativity correctional term, using the positions specified as the parameter.
        Raises an expetion if number_of_bodies!=1
    
    returns:
    arma::Mat<double> positions, a matrix of size dims * number_of_bodies containg the acceleration vectors as column vectors.
    */
    if (number_of_bodies!=1){
        throw("ERROR CAN ONLY ACCEPT ONE BODY");
    }
    arma::Mat<double> forces(dims,1,arma::fill::zeros);
    double l = positions.col(0)(0)*velocities.col(0)(1) - positions.col(0)(1)*velocities.col(0)(0);
    double distA = arma::norm(positions.col(0));
    forces.col(0) = -G1*(positions.col(0))*(1/std::pow(distA,3)+3*l*l/(c*c*std::pow(distA,5)));
    return forces;
};




void SolarSystem::VelocityVerlet(int number_of_timesteps, double dt_length){
    /*
    int number_of_timesteps, the number of timesteps to simulate
    int dt_length, the timstep-length in years

        Integrates the system according to the set up parameters using the VelocityVerlet method.
    */
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
    /*
    int number_of_timesteps, the number of timesteps to simulate
    int dt_length, the timstep-length in years

        Integrates the system according to the set up parameters using the Euler method.
    */
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

void SolarSystem::VelocityVerletMercury(int number_of_timesteps, double dt_length){
     /*
    int number_of_timesteps, the number of timesteps to simulate
    int dt_length, the timstep-length in years

        Integrates the system according to the set up parameters using the Velocity Verlet method using mercury_force as the acceleration.
    */
    dt = dt_length;
    timesteps = number_of_timesteps;
    positions = arma::Cube<double> (dims,number_of_bodies,timesteps,arma::fill::zeros);
    velocities = arma::Cube<double> (dims,number_of_bodies,timesteps,arma::fill::zeros);
    positions.slice(0) = initial_positions;
    velocities.slice(0) = initial_velocities;
    arma::Mat<double> prev_force;
    arma::Mat<double> this_force;
    
    prev_force = mercury_force(positions.slice(0),velocities.slice(0));
    for (int i=1;i<number_of_timesteps;i++){
        positions.slice(i) = positions.slice(i-1) + (velocities.slice(i-1))*dt + (0.5*dt*dt)*prev_force;
        this_force = mercury_force(positions.slice(i),velocities.slice(i-1));
        velocities.slice(i) = velocities.slice(i-1) + (0.5*dt)*(prev_force + this_force);
        prev_force = this_force;
    }
};
void SolarSystem::write_to_file(std::string filename){
    /*
    std::string filename, the file to write the results of the simulations to.

        Write the postions and velocities at every timesteps to filename. Format is specified in the header of the created file.
    */
    std::ofstream output_file;
    output_file.open(filename);

    output_file << "# timesteps: " << timesteps << " dt: " << dt << " nbodies: " << number_of_bodies << " dims: " << dims << "\n";
    output_file << "#each column corresponds to for x11 v11 x21 v21 (x31) (v31) x12 v12 x22 v22 (x32) v(32) ... for each body\n"; 
    output_file << std::setprecision(15);
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
    /*
        Calculates and changes the reference system to the center of mass reference system.
    */
    arma::Col<double> vel_cm = initial_velocities*masses /arma::sum(masses);
    for (int i=0; i<number_of_bodies; i++){
        initial_velocities.col(i) -= vel_cm;
    };
    
};


void SolarSystem::get_init(){
    /*
        
        Prints the initial conditions of the system. Mainly used for debugging purposes.

    */
    std::cout << "Initial Positions:" << std::endl;
    initial_positions.print();
    std::cout << " " << std::endl;
    std::cout << "Initial Velocity:" << std::endl;
    initial_velocities.print();
    std::cout << "\n";
    std::cout << "Masses: " << std::endl;
    masses.print();
};

arma::Cube<double> SolarSystem::get_pos(){
    return positions;
};

arma::Mat<double> SolarSystem::get_init_pos(){
    /*
        
        Return the initial postions of the system. Mainly used for debugging purposes.

    */
    return initial_positions;
};

arma::Mat<double> SolarSystem::get_init_vel(){
    /*
        
        Return the initial velocities of the system. Mainly used for debugging purposes.

    */
    return initial_velocities;
};




SolarSystem::~SolarSystem(){
    /*
    Destructor for the class.
    */

};


/*
//LEGACY CODE: FUNCTION IMPLEMENTATION OF 3A: All the functions here are reused and adapted into the class above. This was only used in development and not in any simulations in the article.
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