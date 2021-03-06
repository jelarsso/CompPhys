#include "ising_model.hpp"

int main(int argc, char* argv[]){
    double start_temp,stop_temp,step_temp;
    int nmc,equiltime;
    std::string filename;

    if (argc != 7){
        std::cout << "Bad usage! This program takes 6 params";
        std::cout << "\n name of output file, number of monte carlo cycles, equiltime, start temperature, stop temperature, and step temperature\n";
        return 1;
    }else{
        filename = argv[1];
        nmc = std::atoi(argv[2]);
        equiltime = std::atoi(argv[3]);
        start_temp = std::stod(argv[4]);
        stop_temp = std::stod(argv[5]);
        step_temp = std::stod(argv[6]);
    }

    IsingModel mdl(20, filename,true,false);
    mdl.Metropolis(nmc,equiltime,start_temp,stop_temp,step_temp);
    return 0;
}