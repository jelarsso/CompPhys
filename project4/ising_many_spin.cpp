#include "ising_model.hpp"

int main(int argc, char* argv[]){
    double start_temp,stop_temp,step_temp;
    int nmc,L;
    std::string filename;

    if (argc != 7){
        std::cout << "Bad usage! This program takes five params";
        std::cout << "\n name of output file, number of monte carlo cycles, L, start temperature, stop temperature, and step temperature\n";
        return 1;
    }else{
        filename = argv[1];
        nmc = std::atoi(argv[2]);
        L = std::atoi(argv[3]);
        start_temp = std::stod(argv[4]);
        stop_temp = std::stod(argv[5]);
        step_temp = std::stod(argv[6]);
    }

    IsingModel mdl(L, filename,false,false);
    mdl.Metropolis(nmc,1000,start_temp,stop_temp,step_temp);
    return 0;
}