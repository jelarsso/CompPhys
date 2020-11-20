#include "ising_model.hpp"

int main(){

    IsingModel mdl(20, "output1.data", false, true);
    for (int N=1; N<1e3+1; N+=1) {
        mdl.Metropolis(N,2.4);
    }
    return 0;
}