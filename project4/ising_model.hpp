#include<armadillo>
#include<string>
#include<fstream>

#ifndef ISING_MODEL_H
#define ISING_MODEL_H


class IsingModel{
    private:
    arma::Mat<int> spin_matrix;
    arma::Col<double> energy_differences;
    double expectation_values[5];
    double initial_temp, final_temp, temp_step, last_temp; 
    double Energy, Magnetization;
    int n_spins;
    int n_mc_cycles,equiltime;
    unsigned long long int accepted_configs;
    bool energy_logger, random_conf;
    std::ofstream output_file;
    std::string filename;
    void Init();


    inline int period(int index, int size);
    void find_energy_differences(double temperature);
    public:
    IsingModel(int number_of_spins, std::string filename, bool el, bool rc);
    ~IsingModel();
    void Metropolis(int number_of_mc_cycles, int etime, double start_temp, double stop_temp, double step_stemp);
    void Metropolis(int number_of_mc_cycles, double temperature);
    void output(double temperature);
};




#endif