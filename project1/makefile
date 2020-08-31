all: compile execute plot

compile:
	c++ -o main.out $(wildcard *.cpp)


execute:
	./main.out

plot:
	python read_data_file.py
