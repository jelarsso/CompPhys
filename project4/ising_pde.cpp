#include "ising_model.hpp"

int main(int argc, char* argv[]){
    double start_temp,stop_temp,step_temp;
    int nmc;
    std::string filename;

    if (argc != 6){
        std::cout << "Bad usage! This program takes 5 params";
        std::cout << "\n name of output file, number of monte carlo cycles, start temperature, stop temperature, and step temperature\n";
        return 1;
    }else{
        filename = argv[1];
        nmc = std::atoi(argv[2]);
        start_temp = std::stod(argv[3]);
        stop_temp = std::stod(argv[4]);
        step_temp = std::stod(argv[5]);
    }

    IsingModel mdl(2, filename,true,false);
    mdl.Metropolis(nmc,start_temp,stop_temp,step_temp);
    return 0;
}