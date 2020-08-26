#include "LUdcmp.hpp"
#include<iostream>
#include<cmath>


using namespace std;


int main(){

    int n = 11;
    double h = 1/double(n+1);

    double a_array[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double b_array[11] = {2,2,2,2,2,2,2,2,2,2,2};
    double c_array[10] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double f_array[11];

    for (int i=0;i<n;i++){
        double x = double(i)/double(n);
        f_array[i] = 100*pow(2.7218281828,-10*x)*h*h;
        cout << f_array[i] << " " ;
    }
    cout << h << endl;

    LUdecomposition lud(n, &a_array[0], &b_array[0], &c_array[0], &f_array[0]);

    lud.findLU();

    lud.print_LU();

    lud.solve();

    lud.print_solution();

    return 0;
}