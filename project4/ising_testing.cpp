#include "ising_model.hpp"

int main(){

    IsingModel mdl(20, "output1.data");
    for (int N=1; N<1e4+1; N+=100) {
        mdl.Metropolis(N,0.1);
    }
    return 0;
}