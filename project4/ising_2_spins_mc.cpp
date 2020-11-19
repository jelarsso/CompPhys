#include "ising_model.hpp"

int main(int argc, char* argv[]){
    double temp;
    int nmc;
    std::string filename;

    if (argc != 4){
        std::cout << "Bad usage! This program takes five params";
        std::cout << "\n name of output file, number of monte carlo cycles, temperature\n";
        return 1;
    }else{
        filename = argv[1];
        nmc = std::atoi(argv[2]);
        temp = std::stod(argv[3]);
    }

    IsingModel mdl(2, filename);
    mdl.Metropolis(nmc,temp);
    return 0;
}