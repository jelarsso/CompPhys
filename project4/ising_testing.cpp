#include "ising_model.hpp"

int main(){
    IsingModel mdl(20, "output8.data", false, true);
    for (int N=1; N<1e4+1; N+=10) {
        mdl.Metropolis(N,2.4);
    }
    return 0;
}