#include<string>
#include<armadillo>
#ifndef SOLAR_SYSTEM_H
#define SOLAR_SYSTEM_H


class SolarSystem{
    private:
    arma::Cube<double>* positions;
    arma::Cube<double>* velocities;
    arma::Mat<double>* initial_positions;
    arma::Mat<double>* initial_velocities;
    arma::Row<double>* masses;
    int number_of_bodies;
    int dims;

    void Verlet_step(int number_of_timesteps, int dt);
    void Euler_step(int number_of_timesteps, int dt);
    arma::Mat<double> force(arma::Mat<double> positions);

    public:
    SolarSystem(int number_of_bodies, int dims);
    ~SolarSystem();
    void read_state(std::string filename);
    void solve(int number_of_timesteps, int dt, bool method);
    void write_result(std::string filename);
};



#endif