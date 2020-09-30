# Description of code for project 2


## Jacobi rotation algorithm
The jacobi rotation is implemented in jacobi_rotation.cpp and its header jacobi_rotation.hpp.
The functions are documented in jacobi_rotation.cpp.

## Running the code
The c++ code for the different parts of the project was mainly executed from python. We made all plots and timed the programs from python. The python files only requires that the c++ programs have been compiled with the specified executable names given in the makefile.
The file compare_solution.py runs and plots the results needed for part 2b.
The file run_project.py implements several functions which runs and solves a specific experiment, see the comments on the functions for more details.

## Makefile
The makefile can compile all the required programs, however the path to armadillo must be specified for your computer.
Run with:

```    make ```

## Buckling beam / fixed spring
The problem corresponding to the buckling beam or fixed spring is implemented in main.cpp. It prints the eigenvalues and writes them and the eigenvectors to a file.
Compile the code with:

```    c++ main.cpp jacobi_rotation.cpp -o main.o ```

Run the code with

```    ./main.o size tol ```

where size is the an integer size of the matrix and tol the tolerance for the jacobi algorithm.

## Quantum dots one/two electron
The two files quantum_dots_(one/two)_electron.cpp corresponds to the solution for the one and two electron quantum mechanical systems described in the article.
Compile the code with:

```    c++ quantum_dots_(one/two)_electron.cpp jacobi_rotation.cpp -o quant.o ```

Run the code with:

```    ./quant.o size tol omega_r ```

where size is the integer size of the matrix, tol the tolerance and (only for the two electron system) omega_r the strength of the harmonic potential.

## Armadillo implementation
The armadillo implementation in armadillo_implementatio.cpp solves the same problem as main.cpp.
Compile with:

```    g++ armadillo_implementation.cpp -o arma.o -O2 -I /path/to/armadillo/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack ```

Run the code with:

```    ./arma.o size ```

where size is the an integer size of the matrix.

## Testing
Unit tests have been implemented in test_jacobi_rotation.cpp and is compiled with:

 ```   g++ test_jacobi_rotation.cpp jacobi_rotation.cpp -o test.o -O2 -I /path/to/armadillo/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack```

make sure that the catch2 (from https://github.com/catchorg/Catch2) header file catch.hpp is in the same directory as the test_jacobi_rotation.cpp file.
The test are then run with:

```    ./test.o ```

Or if the makefile has been run:

```    make test ```

Our tested result yields,

```\>>>.-/test.o
===============================================================================
All tests passed (10101 assertions in 1 test case)```

