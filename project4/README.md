# Description of Code for Project 4


# Timing and optimalization of ising_parallel.cpp
We parallelized the code using openmp.
We mainly parallelized the outer temperature loop, that is independent MC simulations for different temperatures run concurrently, thus the biggest benefit from parallelization is only seen when simulating many different temperatures.
The following results show the performace increase for a select simulation.

Compiled using
c++ ising_parallel.cpp -o para -fopenmp -O3  -I /home/Documents/armadillo-9.900.3/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

Ran on AMD Ryzen 5 3500U with 8 threads 3.7 GHz

Without paralellization
(base) johan@johanemil:~/Documents/UiO/CompPhys/project4$ time ./para "output.data" 100000 1000 20 2.0 2.6 0.05
real	0m25,206s

With parallelization
(base) johan@johanemil:~/Documents/UiO/CompPhys/project4$ time ./para "output.data" 100000 1000 20 2.0 2.6 0.05
real	0m5,048s
