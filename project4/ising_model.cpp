#include "ising_model.hpp"
#include <armadillo>
#include <random>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

inline int IsingModel::period(int index, int size){
    /*
    In:
    int index, index to index the array 
    int size, size of the array
    
    Finds the correct index to an array with periodic boundary conditions.
    
    Out:
    int the new index
    */
    return (index+size)%size;
};

IsingModel::IsingModel(int number_of_spins, std::string filen, bool el, bool rc){
    /*
    In:
    int number_of_spins, the grid size or number of spins L*L'
    string filename, the filename
    bool el, rc, whether to turn on or off the energy logging or random start configurations
    
    Constructor.
    */
    n_spins = number_of_spins;
    filename = filen;
    energy_logger = el;
    random_conf = rc;
};

void IsingModel::Init(){
    /*
    Initialize the array of the spins according to the parameters set by the constructor.
    */
    if (!random_conf){
    spin_matrix.ones(n_spins,n_spins);
    }else{
        std::random_device rd1;
        std::mt19937 l(rd1());
        std::uniform_int_distribution<int> distribution(0,1);

        spin_matrix.ones(n_spins,n_spins);
        for (int x = 0; x<n_spins; x++){
            for (int y = 0; y<n_spins; y++){
                int a = distribution(l);
                spin_matrix(x,y) = a*2-1;
            };
        };
    }
    Magnetization = (double) arma::accu(spin_matrix);
    Energy = 0;
    accepted_configs = 0;
    for (int x = 0; x<n_spins; x++){
        for (int y = 0; y<n_spins; y++){
            Energy -= (double) spin_matrix(x,y)*(spin_matrix(period(x-1,n_spins),y)+spin_matrix(x,period(y-1,n_spins)));
        }
    }
}


IsingModel::~IsingModel(){
    /*
    Destructor for the class. 
    */
    output_file.close();
};

void IsingModel::find_energy_differences(double temperature){
    /*
    double temperature: the temperature to calculate the differences at.
    Finds the energy differences for the possible spin transistions at temperature.
    */
    energy_differences.zeros(17);
    for (int deltaE = -8; deltaE<=8; deltaE+=4){
        energy_differences(deltaE+8) = std::exp(-deltaE/temperature);
    }  
};

void IsingModel::Metropolis(int number_of_mc_cycles, int etime, double start_temp, double stop_temp, double temperature_step){
    /*
    int number_of_mc_cycles: The number of monte carlo cycles to perform.
    int etime: The equilibration time. The number of MC cycles to perform before starting to calculate the expectation values.
    double start_temp: The temperature to start the simulations at.
    double stop_temp: The temperature to end the simulations at.
    double temperature_step: The step between the temperatures between start and stop.
    
    Simulate the Ising model with the specified parameters. Will perform several simulations between start_temp and stop_temp.
    Outputs to filename
    */
    n_mc_cycles = number_of_mc_cycles;
    initial_temp = start_temp;
    final_temp = stop_temp;
    temp_step = temperature_step;
    equiltime = etime;
    for (double temp = initial_temp; temp<=final_temp; temp+=temp_step){
        Metropolis(number_of_mc_cycles,temp);
    }
};


void IsingModel::Metropolis(int number_of_mc_cycles, double temp){
    /*
    int number_of_mc_cycles: the number of MC cycles
    double temp: the temperature to simulate at.
    
    Perform the Metropolis algorithm for one temperature. 
    */
    n_mc_cycles = number_of_mc_cycles;
    last_temp = temp;

    //Seed Mersenne twister
    std::random_device rd;
    std::mt19937_64 gen(rd());

    std::uniform_real_distribution<double> distribution(0.0,1.0);

    Init();
    find_energy_differences(temp);
    for (int i=0;i<5;i++) expectation_values[i]=0;

    if (energy_logger==true){
        if (output_file.is_open() == false){
        output_file.open(filename);
        output_file << "#n_spins " << n_spins << " n_mc_cycles " << n_mc_cycles << " all values per spin\n";
        output_file << "# temp Eavg  Evar  Mavg  Mvar  Mabsavg accepted_configs\n"; 
        }
    output_file << "# Energy logger is ON, performance is reduced, the next line is the energy at every monte carlo cycle \n";
    }

    for (int cycle = 0; cycle<n_mc_cycles; cycle++){
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
    if (energy_logger==true && cycle>equiltime){
    output_file << Energy/n_spins/n_spins << " ";
    }
    if (cycle>equiltime){
    expectation_values[0] += Energy;
    expectation_values[1] += Energy*Energy;
    expectation_values[2] += Magnetization;
    expectation_values[3] += Magnetization*Magnetization;
    expectation_values[4] += std::fabs(Magnetization);
    }
    }
    if (energy_logger==true){
       output_file << "\n";
    }
    output(temp);
};



void IsingModel::output(double temperature){
    /*
    double temperature: the temperature that the last simulation was performed at, to get correct normalization.
    
    Function to output the expectation values to filename.
    */
    double norm = 1/(double) (n_mc_cycles - equiltime);
    double Eaverage = expectation_values[0]*norm;
    double E2average = expectation_values[1]*norm;
    double Maverage = expectation_values[2]*norm;
    double M2average = expectation_values[3]*norm;
    double Mabsaverage = expectation_values[4]*norm;
    double Evariance = (E2average - Eaverage*Eaverage)/n_spins/n_spins;
    double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
    
    if (output_file.is_open() == false){
        output_file.open(filename);
        output_file << "#n_spins " << n_spins << " n_mc_cycles " << n_mc_cycles << " equiltime " << equiltime << " all values per spin\n";
        output_file << "# temp Eavg  Evar  Mavg  Mvar  Mabsavg accepted_configs\n"; 
    }
    
    output_file << std::setprecision(15);
    output_file << temperature << " ";
    output_file << Eaverage/n_spins/n_spins << " ";
    output_file << Evariance/temperature/temperature << " ";
    output_file << Maverage/n_spins/n_spins << " ";
    output_file << Mvariance/temperature << " ";
    output_file << Mabsaverage/n_spins/n_spins << " ";
    output_file << accepted_configs << "\n";
};