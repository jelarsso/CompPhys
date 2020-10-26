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
