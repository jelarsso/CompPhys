#include "ising_model.hpp"

int main(){

    IsingModel mdl(2, "output1.data");
    mdl.Metropolis(1000,1,2,0.1);
    return 0;
}