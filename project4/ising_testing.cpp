#include "ising_model.hpp"

int main(){

    IsingModel mdl(20, "output1.data");
    for (int N=1; N<101; N++) {
        mdl.Metropolis(N,1);
    }
    return 0;
}