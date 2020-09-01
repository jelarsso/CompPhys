all: compile execute plot

compile:
	c++ -o main.out $(wildcard *.cpp)

special:
	c++ -o special.out specialized_algorithm.cpp LUdcmp.cpp

general:
	c++ -o general.out general_algorithm.cpp LUdcmp.cpp


execute:
	./main.out

plot:
	python read_data_file.py
