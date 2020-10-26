#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "solar_system.hpp"
#include<armadillo>
#include<cmath>

const double pi = 3.14159265358979;


TEST_CASE() {
    //Test that the intial conditions has the correct dimensionalities and values.
    int *indices;
    indices = new int[2];
    indices[0] = 3;
    indices[1] = 5;

    double simulation_length = 3;
    double dt = 1e-4;
    double start_vel = 2*pi;
    double beta = 2; 
    int timesteps = (int)(simulation_length/dt)+1;
    int *indices;
    indices = new int[1];
    indices[0] = 1;

    arma::Mat<double> init_pos = arma::Mat<double>(2,1);
    arma::Mat<double> init_vel = arma::Mat<double>(2,1);
    arma::Col<double> masses = arma::Col<double>(1);
    init_pos(0,0) = 1;
    init_pos(1,0) = 0;
    init_vel(0,0) = 0;
    init_vel(1,0) = start_vel;
    masses(0) = 1;
    

    SolarSystem sol(2,1,indices,beta,false,true);
    sol.set_initial_conditions(init_pos,init_vel,masses);
    sol.VelocityVerlet(timesteps,dt);

    arma::Cube<double> position = sol.get_pos();
    for (int i = 0; i<timesteps; i++){
        REQUIRE(std::abs(position(0, 0, i)*position(0, 0, i) + position(1, 0, i)*position(1, 0, i)) - 1 < 0.01);
    }
    
    SolarSystem sol(3,2,indices,2,false,true);
    arma::Mat<double> pos;
    pos = sol.get_init_pos();
    arma::Mat<double> vel;
    vel = sol.get_init_vel();
    sol.get_init();
    pos.print();

    REQUIRE(pos.n_cols==2);
    REQUIRE(pos.n_rows==3);
    REQUIRE(vel.n_cols==2);
    REQUIRE(vel.n_rows==3);

    REQUIRE(pos(0,0)==9.323412692535679E-01);
    REQUIRE(vel(2,1)==-1.616464625977914E-04*365);
};
