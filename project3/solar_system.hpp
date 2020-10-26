#include<string>
#include<armadillo>
#ifndef SOLAR_SYSTEM_H
#define SOLAR_SYSTEM_H


class SolarSystem{
    private:
    int dims;
    int timesteps;
    double simulation_length;
    double dt;
    double beta;
    int number_of_bodies;
    int* body_indices;
    bool stationary_sun;

    arma::Mat<double> initial_positions;
    arma::Mat<double> initial_velocities;
    arma::Col<double> masses;

    arma::Cube<double> positions;
    arma::Cube<double> velocities;
    
    void read_initial_condition(std::string filename,int body_index, double* mass,arma::Col<double>* initial_position, arma::Col<double>* initial_velocity);
    arma::Mat<double> force(arma::Mat<double> positions);
    arma::Mat<double> stat_sun_force(arma::Mat<double> positions);
    arma::Mat<double> mercury_force(arma::Mat<double> positions, arma::Mat<double> velocities);
    void read_initial_conditions(std::string filename);

    public:
    SolarSystem(int dims, int number_of_bodies, int* body_indices, double beta, bool set_initial_conditions_manually, bool stationary_sun);
    ~SolarSystem();

    void get_init();
    arma::Mat<double> get_init_pos();
    arma::Mat<double> get_init_vel();
    void set_initial_conditions(arma::Mat<double> initial_positions, arma::Mat<double> initial_velocities, arma::Col<double> mass);
    void set_initial_conditions(arma::Col<double> mass);

    arma::Cube<double> get_pos();

    void Euler(int number_of_timesteps, double dt);
    void VelocityVerlet(int number_of_timesteps, double dt);
    void VelocityVerletMercury(int number_of_timesteps, double dt_length);


    void write_to_file(std::string filename);
    
    void change_reference_to_cm();
    
};

/*LEGACY CODE: Function implementations.
void read_inital_condition(std::string filename,int body_index, double* mass, arma::Col<double>* initial_position, arma::Col<double>* initial_velocity);
void read_initial_conditions(std::string filename, int number_of_bodies, int body_indices[], arma::Col<double>* masses, arma::Mat<double>* initial_position, arma::Mat<double>* initial_velocity);

void Euler(int number_of_timesteps, double dt, int number_of_bodies, arma::Col<double> masses, arma::Mat<double>* initial_positions, arma::Mat<double>* initial_velocities, arma::Cube<double>* positions, arma::Cube<double>* velocities);
void VelocityVerlet(int number_of_timesteps, double dt, int number_of_bodies, arma::Col<double> masses, arma::Mat<double>* initial_positions, arma::Mat<double>* initial_velocities, arma::Cube<double>* positions, arma::Cube<double>* velocities);

arma::Mat<double> force(int number_of_bodies,arma::Col<double> masses, arma::Mat<double> positions);

void write_to_file(int dims, int number_of_bodies, int number_of_timesteps, arma::Cube<double> positions, std::string filename);
*/

#endif