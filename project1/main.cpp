#include "LUdcmp.hpp"
#include<iostream>
#include<cmath>


using namespace std;


int main(){

    int n = 11;
    double h = 1/double(n-1);

    double f_array[n];

    for (int i=0;i<n;i++){
        double x = double(i)*h;
        f_array[i] = 100*exp(-10*x)*h*h;
        }

    LUdecomposition lud(n, &f_array[0]);

    lud.findLU();

    lud.print_LU();

    lud.solve();

    lud.print_solution();

    return 0;
}