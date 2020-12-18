# Project 5: The Diffusion Project

We chose the Diffusion Project and did the optional Lithosphere part.
This folder contains the code for the project in divided in the following files.
Make sure to edit the makefile to point the gcc-compiler to the correct location for the armadillo-library.

Requirements are recent versions of python (>3.7) and the module plotly along with the scipy-stack.
And a c++-compiler with armadillo (>9.900.2)

algorithms.cpp - contains the functions implemented in .cpp and its .hpp header file.

p5c.cpp - Uses the 3 algorithms to solve the one dimensional problem described in the report / task. It also calculates the analytical result.
    It is compiled with the command "make 5c" which yields the executable "p5c". To run the program type "./p5c dx nt" where dx is the step length and nt the number of time points. It outputs 4 files "analytical.data", "backward_euler.data", "forward_euler.data" and "cnicholson.data". The program also prints nx, the number of grid points, such that reading the file and using np.loadtxt(file).reshape((nt+1, nx+1)) in python loads the simulated values or analytical values into an easily accessible np-array. 

p5f.cpp - This uses the Forward Euler algorithm in two dimensions.
    Compile it using "make 5f". Run with "./p5f dx nt" as above, which yields a nx+1 by nx+1 grid solves it and writes it to the data file "fe2d.data" along with the analytical result in "analytical2d.data". To access in python simply use the following line of code: "np.loadtxt(file).reshape((nt+1,nx+1,nx+1))".


litho.cpp - This simulates the lithosphere in 3 parts.
    First the non-radioactive state, then the pre-enriched and the enriched, and writes it to threes files "litho_no_Q.data", "litho_no_Q_pb.data" and "litho_enriched.data". All scaling of variables is done in this file. Compile "make lithos", run "./litho dx T" where dx is the step length in x' units as described in the text and T the number of Gy to simulate for.


main.py - The python file that does the plotting and analysis of the data files from c++. See the file for more details. If numpy fails to resize the array, make sure that the resize-command has the corresponding nt,nx,ny values that the c++-program that had when executed.