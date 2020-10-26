# Description of code.

The main files are the implementation of the Solar system class in solar_system.(c/h)pp and the python-file data_analyzer.py which is used to execute the scripts and generate the plots.

A description of the class is given in the solar_system.cpp-file.

The compiled programs used in the python script are the following files:
~ earth_sun_system.cpp
~ mercury_sim.cpp
~ three_body.cpp
They are each described in their file.

The programs are compiled on linux with the simple command,
~ c++ earth_sun_system.cpp solar_system.cpp -o earth_sun_system -O2 -I /path/to/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack

The test in test_solar_system.cpp are compiled with,
~ c++ test_solar_system.cpp solar_system.cpp catch.hpp -o test -O2 -I /path/to/armadillo-9.900.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack