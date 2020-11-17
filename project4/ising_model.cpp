#include "ising_model.hpp"
#include <armadillo>
#include <random>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

int IsingModel::period(int index, int size){
    return (index+size)%size;
};

IsingModel::IsingModel(int number_of_spins, std::string filen){
    n_spins = number_of_spins;
    filename = filen;
    Init();
};

void IsingModel::Init(){
    spin_matrix.ones(n_spins,n_spins);
    Magnetization = (double) arma::accu(spin_matrix);
    Energy = 0;
    for (int x = 0; x<n_spins; x++){
        for (int y = 0; y<n_spins; y++){
            Energy += (double) spin_matrix(x,y)*(spin_matrix(period(x-1,n_spins),y)*spin_matrix(x,period(y-1,n_spins)));
        }
    }
}


IsingModel::~IsingModel(){
    output_file.close();
};

void IsingModel::find_energy_differences(double temperature){
    energy_differences.zeros(17);
    for (int deltaE = -8; deltaE<=8; deltaE+=4){
        energy_differences(deltaE+8) = std::exp(-deltaE/temperature);
    }  
};

void IsingModel::Metropolis(int number_of_mc_cycles, double start_temp, double stop_temp, double temperature_step){
    n_mc_cycles = number_of_mc_cycles;
    initial_temp = start_temp;
    final_temp = stop_temp;
    temp_step = temperature_step;
    for (double temp = initial_temp; temp<=final_temp; temp+=temp_step){
        Metropolis(number_of_mc_cycles,temp);
    }
};


void IsingModel::Metropolis(int number_of_mc_cycles, double temp){
    n_mc_cycles = number_of_mc_cycles;
    last_temp = temp;

    //Seed Mersenne twister
    std::random_device rd;
    std::mt19937_64 gen(rd());

    std::uniform_real_distribution<double> distribution(0.0,1.0);

    Init();
    find_energy_differences(temp);
    expectation_values.zeros(5);

    for (int cycle = 0; cycle<n_mc_cycles; cycle++){
        for (int x = 0; x<n_spins; x++){
            for (int y = 0; y<n_spins; y++){
                int ix = (int) (distribution(gen)*(double)n_spins);
                int iy = (int) (distribution(gen)*(double)n_spins);
                int deltaE = -2*spin_matrix(ix,iy)*(spin_matrix(period(ix+1,n_spins),iy)+spin_matrix(period(ix-1,n_spins),iy)+spin_matrix(ix,period(iy+1,n_spins))+spin_matrix(ix,period(iy-1,n_spins)));

                if (distribution(gen) <= energy_differences(deltaE+8)){
                    spin_matrix(ix,iy)*=-1;
                    Magnetization += (double) 2*spin_matrix(ix,iy);
                    Energy += (double) deltaE;
                }
            }
        }
    expectation_values(0) += Energy;
    expectation_values(1) += Energy*Energy;
    expectation_values(2) += Magnetization;
    expectation_values(3) += Magnetization*Magnetization;
    expectation_values(4) += std::abs(Magnetization);
    }

    output(temp);
};



void IsingModel::output(double temperature){
    double norm = 1/(double) n_mc_cycles;
    double Eaverage = expectation_values(0)*norm;
    double E2average = expectation_values(1)*norm;
    double Maverage = expectation_values(2)*norm;
    double M2average = expectation_values(3)*norm;
    double Mabsaverage = expectation_values(4)*norm;
    double Evariance = (E2average- Eaverage*Eaverage);
    double Mvariance = (M2average - Mabsaverage*Mabsaverage);
    
    if (output_file.is_open() == false){
        output_file.open(filename);
        output_file << "#n_spins" << n_spins << " n_mc_cycles " << n_mc_cycles << "\n";
        output_file << "# temp Eavg  Evar  Mavg  Mvar  Mabsavg\n"; 
    }

    output_file << std::setprecision(15);
    output_file << temperature << " ";
    output_file << Eaverage << " ";
    output_file << Evariance/temperature/temperature << " ";
    output_file << Maverage << " ";
    output_file << Mvariance/temperature << " ";
    output_file << Mabsaverage << "\n";
};