#include "LUdcmp.hpp"
#include<iostream>
#include<cmath>
#include<string>

using namespace std;

int main(){

    int n = 10000;
    double h = 1/double(n-1);

    //double f_array[n];
    double *f_array;
    f_array = new double[n];

    for (int i=0;i<n;i++){
        double x = double(i)*h;
        f_array[i] = 100*exp(-10*x)*h*h;
        }


    LUdecomposition lud(n, &f_array[0]);

    //lud.findLU();

    //lud.print_LU();

    lud.solve();

    string filename = "output.data";
    lud.write_v_to_file(filename);

    //lud.print_solution();

    return 0;
}