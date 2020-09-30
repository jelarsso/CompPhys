compile:
	c++ main.cpp jacobi_rotation.cpp -o main.o
	c++ quantum_dots_one_electron.cpp jacobi_rotation.cpp -o quant.o
	c++ quantum_dots_two_electron.cpp jacobi_rotation.cpp -o quant_two.o
	c++ armadillo_implementation.cpp -o arma.o -O2 -I /home/johan/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack
	c++ test_jacobi_rotation.cpp jacobi_rotation.cpp -o test.o -O2 -I /home/johan/armadillo/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

test:
	./test.o