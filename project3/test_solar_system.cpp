#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "solar_system.hpp"
#include<armadillo>
#include<cmath>

TEST_CASE() {
    //Test that the intial conditions has the correct dimensionalities and values.
    int *indices;
    indices = new int[2];
    indices[0] = 3;
    indices[1] = 5;

    SolarSystem sol(3,2,indices,2,false,true);
    arma::Mat<double> pos = sol.get_init_pos();
    arma::Mat<double> vel = sol.get_init_vel();
            
};
