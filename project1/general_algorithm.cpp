#include "LUdcmp.hpp"
#include<iostream>
#include<cmath>
#include<string>

using namespace std;

int main(int argc, char* argv[]){
    int n = atoi(argv[1]);
    double h = 1/double(n);

    //solve the system using the general algorithm

    double *f_array;
    double *a_array;
    double *b_array;
    double *c_array;
    f_array = new double[n];
    a_array = new double[n-1];
    b_array = new double[n];
    c_array = new double[n-1];
    

    for (int i=0;i<n;i++){
        double x = double(i)*h;
        f_array[i] = 100*exp(-10*x)*h*h;
        b_array[i] = 2;
        }
    
    for (int i=0;i<n-1;i++){
        a_array[i] = -1;
        c_array[i] = -1;
    }


    LUdecomposition lud(n, &a_array[0], &b_array[0], &c_array[0],  &f_array[0]);

    lud.findLU();
    lud.solve();
    string filename = "general_output.data";
    lud.write_v_to_file(filename);
    delete[] f_array;
    delete[] a_array;
    delete[] b_array;
    delete[] c_array;
    return 0;
}