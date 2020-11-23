# Description of Code for Project 4

The main program is the main.py python script which calls the required C++ programs and reads and plots the output.

Run with:
> python main.py

To compile the C++, run 
> make project

The makefile may need to be adapted to add the correct path to the armadillo library.

For example if you wish to run the Ising model for a 50*50 grid, the following call could be used.
> ./para filename nmc equiltime L start_temp stop_temp step_temp

where the filename is the filename, nmc is the number of monte carlo cycles, equiltime the equilibration time (note it must be less than nmc), L denotes the grid size (L*L), starting at start_temp and ending at stop_temp with step_temp step length.


# Timing and optimalization of ising_parallel.cpp
We parallelized the code using openmp.
We mainly parallelized the outer temperature loop, that is independent MC simulations for different temperatures run concurrently, thus the biggest benefit from parallelization is only seen when simulating many different temperatures. One can note that the initial ising_model.cpp was refactored away from a class and to functions to make the parallellization easier.
The following results show the performace increase for a select simulation.

Compiled using
> c++ ising_parallel.cpp -o para -fopenmp -O3  -I /home/Documents/armadillo-9.900.3/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

Ran on AMD Ryzen 5 3500U with 8 threads 3.7 GHz

Without paralellization
> (base) johan@johanemil:~/Documents/UiO/CompPhys/project4$ time ./para "output.data" 100000 1000 20 2.0 2.6 0.05

> real	0m25,206s

With parallelization
> (base) johan@johanemil:~/Documents/UiO/CompPhys/project4$ time ./para "output.data" 100000 1000 20 2.0 2.6 0.05

> real	0m5,048s
