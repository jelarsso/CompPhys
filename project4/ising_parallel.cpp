#include <armadillo>
#include <random>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

//Definitions
std::ofstream output_file;
inline int period(int index, int size){
    return (index+size)%size;
};
void Init(arma::Mat<int>*, double*, double*, int*, int);
arma::Col<double> find_energy_differences(double);
void output(double*, int, std::string, double, int, int, int);
void Metropolis(std::string, int, int, int, double);
void IsingModelMetropolis(std::string, int, int, int, double, double, double);

//Implementations

/*
The functions in this file is copied from ising_model.cpp but the class structure has been removed.
The int main() function shows the usage of the function, however it is simply needed to compile with "make parallel"
*/

void Init(arma::Mat<int>* spin_matrix, double* Energy, double* Magnetization, int* accepted_configs, int n_spins){
    spin_matrix->ones(n_spins,n_spins);
    *Magnetization = (double) arma::accu(*spin_matrix);
    *Energy = 0;
    *accepted_configs = 0;
    for (int x = 0; x<n_spins; x++){
        for (int y = 0; y<n_spins; y++){
            *Energy -= (double) (*spin_matrix)(x,y)*((*spin_matrix)(period(x-1,n_spins),y)+(*spin_matrix)(x,period(y-1,n_spins)));
        }
    }
}



arma::Col<double> find_energy_differences(double temperature){
    arma::Col<double> energy_differences(17,arma::fill::zeros);
    for (int deltaE = -8; deltaE<=8; deltaE+=4){
        energy_differences(deltaE+8) = std::exp(-deltaE/temperature);
    }
    return energy_differences;
};


void output(double* expectation_values, int accepted_configs, std::string filename, double temperature, int nmccycles, int equiltime, int n_spins){
    double norm = 1/((double) (nmccycles-equiltime));
    double Eaverage = expectation_values[0]*norm;
    double E2average = expectation_values[1]*norm;
    double Maverage = expectation_values[2]*norm;
    double M2average = expectation_values[3]*norm;
    double Mabsaverage = expectation_values[4]*norm;
    double Evariance = (E2average - Eaverage*Eaverage)/n_spins/n_spins;
    double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
    
    #pragma omp critical
    {
    output_file << std::setprecision(15);
    output_file << temperature << " ";
    output_file << Eaverage/n_spins/n_spins << " ";
    output_file << Evariance/temperature/temperature << " ";
    output_file << Maverage/n_spins/n_spins << " ";
    output_file << Mvariance/temperature << " ";
    output_file << Mabsaverage/n_spins/n_spins << " ";
    output_file << accepted_configs << "\n";
    }
};


void Metropolis(std::string filename, int n_spins, int nmccycles, int equiltime, double temp){

    //Seed Mersenne twister
    std::random_device rd;
    std::mt19937_64 gen(rd());

    std::uniform_real_distribution<double> distribution(0.0,1.0);

    double Energy, Magnetization;
    int accepted_configs;
    arma::Mat<int> spin_matrix;

    Init(&spin_matrix,&Energy,&Magnetization,&accepted_configs,n_spins);
    arma::Col<double> energy_differences = find_energy_differences(temp);
    double expectation_values[5];
    for (int i=0;i<5;i++) expectation_values[i]=0;

    
    for (int cycle = 0; cycle<nmccycles; cycle++){
        for (int x = 0; x<n_spins; x++){
            for (int y = 0; y<n_spins; y++){
                int ix = (int) (distribution(gen)*(double)n_spins);
                int iy = (int) (distribution(gen)*(double)n_spins);
                int deltaE = 2*spin_matrix(ix,iy)*(spin_matrix(period(ix+1,n_spins),iy)
                +spin_matrix(period(ix-1,n_spins),iy)
                +spin_matrix(ix,period(iy+1,n_spins))
                +spin_matrix(ix,period(iy-1,n_spins)));

                if (distribution(gen) <= energy_differences(deltaE+8)){
                    spin_matrix(ix,iy)*=-1;
                    Magnetization += (double) 2*spin_matrix(ix,iy);
                    Energy += (double) deltaE;
                    accepted_configs++;
                }

            }
        }
    if (cycle>equiltime){
    expectation_values[0] += Energy;
    expectation_values[1] += Energy*Energy;
    expectation_values[2] += Magnetization;
    expectation_values[3] += Magnetization*Magnetization;
    expectation_values[4] += std::fabs(Magnetization);
    }
    }
    output(expectation_values, accepted_configs, filename, temp, nmccycles, equiltime, n_spins);
};

void IsingModelMetropolis(std::string filename, int n_spins, int nmccycles, int equiltime, double start_temp, double stop_temp, double temp_step){
    if (output_file.is_open() == false){
        output_file.open(filename);
        output_file << "#n_spins " << n_spins << " nmccycles " << nmccycles << " all values per spin\n";
        output_file << "# temp Eavg  Evar  Mavg  Mvar  Mabsavg accepted_configs\n"; 
    }
    int N = (stop_temp-start_temp)/temp_step;
    omp_set_num_threads(8);
    #pragma omp parallel for
    for (int i = 0; i<=N; i++){
        double temp = start_temp + i*temp_step;
        Metropolis(filename, n_spins, nmccycles, equiltime, temp);
    }
    output_file.close();
};


int main(int argc, char* argv[]){
    double start_temp,stop_temp,step_temp;
    int nmc,L,equiltime;
    std::string filename;

    if (argc != 8){
        std::cout << "Bad usage! This program takes five params";
        std::cout << "\n name of output file, number of monte carlo cycles, equilibration time, L, start temperature, stop temperature, and step temperature\n";
        return 1;
    }else{
        filename = argv[1];
        nmc = std::atoi(argv[2]);
        equiltime = std::atoi(argv[3]);
        L = std::atoi(argv[4]);
        start_temp = std::stod(argv[5]);
        stop_temp = std::stod(argv[6]);
        step_temp = std::stod(argv[7]);
    }
    IsingModelMetropolis(filename, L, nmc, equiltime, start_temp, stop_temp, step_temp);
    return 0;
}