
2spin:
	c++ ising_2_spins.cpp ising_model.cpp -o spins2 -O3 -I /home/johan/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

nspin:
	c++ ising_many_spin.cpp ising_model.cpp -o spins -O3 -I /home/johan/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

2spinnmc:
	c++ ising_2_spins_mc.cpp ising_model.cpp -o spins2_mc -O3 -I /home/johan/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

pde:
	c++ ising_pde.cpp ising_model.cpp -o pde -O3 -I /home/johan/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

testing:
	c++ ising_testing.cpp ising_model.cpp -o testing -O2 -I /home/johan/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

op:
	c++ ising_testing.cpp ising_model.cpp -o testing -O3 -I /home/Documents/armadillo-9.900.3/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

parallel:
	c++ ising_parallel.cpp -o para -fopenmp -O3  -I /home/Documents/armadillo-9.900.3/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

project:
	2spin nspin pde op parallel