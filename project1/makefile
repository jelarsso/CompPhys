all: compile execute plot

compile:
	c++ -o main.out $(wildcard *.cpp)

special:
	c++ -o special.out specialized_algorithm.cpp LUdcmp.cpp

general:
	c++ -o general.out general_algorithm.cpp LUdcmp.cpp

arma:
	c++ armadillo_algorithm.cpp -o arma.out -O2 -I /home/johan/armadillo/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

execute:
	./main.out

project: general special arma
	python run_code.py
