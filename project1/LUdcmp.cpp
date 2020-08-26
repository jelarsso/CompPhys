#include "LUdcmp.hpp"
#include<iostream>


using namespace std;



LUdecomposition::LUdecomposition(int N, double *f_array){
    n = N;
    h_stepsize = 1/double(n+1);

    constant_abc = true;
    
    a = new double;
    b = new double;
    c = new double;

    *a = -1;
    *b = 2;
    *c = -1;
    
    f = f_array;
};


LUdecomposition::LUdecomposition(int N, double *a_array, double *b_array, double *c_array, double *f_array){
    n = N;
    h_stepsize = 1/double(n+1);

    constant_abc = false;
    
    a = a_array;

    b = b_array;

    c = c_array;

    f = f_array;
};


LUdecomposition::~LUdecomposition(){
    if (constant_abc){
        delete a;
        delete b;
        delete c;
    }
    if (LU_found){
    delete[] d;
    delete[] l;
    }
    if (solved){
    delete[] u;
    delete[] v;
    }
};

double *LUdecomposition::geta(){
    return a;
};

void LUdecomposition::findLU(){
    
    if (constant_abc){
        findLU_constant_abc();
    }else{
        findLU_varying_abc();
    }
    LU_found = true;
};


void LUdecomposition::findLU_constant_abc(){ //not strictly needed
    l = new double[n-1];
    d = new double[n];

    d[0] = *b;

    for (int i = 1; i<n; i++){
        l[i-1] = *a/d[i-1];
        d[i] = *b - l[i-1]*(*c);
    }
};

void LUdecomposition::findLU_varying_abc(){
    l = new double[n-1];
    d = new double[n];

    d[0] = b[0];

    for (int i = 1; i<n; i++){
        l[i-1] = a[i-1]/d[i-1];
        d[i] = b[i] - l[i-1]*c[i-1];
    }
};


void LUdecomposition::print_LU(){
    if (!LU_found){
        throw "findLU not yet called!";
    }

    cout << "d = [";
    for (int i=0;i<n;i++){
        cout << " " << d[i];
    }
    cout << " ]" << endl;


    cout << "l = [";
    for (int i=0;i<n-1;i++){
        cout << " " << l[i];
    }
    cout << " ]" << endl;
};

void LUdecomposition::solve(){
    
    if (constant_abc){
        solve_constant_abc();
    }else{
        solve_varying_abc();
    }
    solved = true;
};


void LUdecomposition::solve_constant_abc(){
    u = new double[n];
    v = new double[n];

    u[0] = *b;

    for (int i = 1; i<n; i++){
        u[i] = 2 + u[i-1]*double(i-1)/double(i);
    }
    v[n] = -double(n+1)/double(n);
    for (int i=n-1; i>0; i--){
        v[i] = (u[i] + v[i+1])*float(i)/float(i+1);
    }
};

void LUdecomposition::solve_varying_abc(){
    if (!LU_found){
        throw "LU-decomposition must be found before calling this function";
    }
    u = new double[n];
    v = new double[n];

    u[0] = b[0];

    for (int i = 1; i<n; i++){
        u[i] = b[i] - u[i-1]*l[i-1];
    }
    v[n]  = d[n]/a[n];

    for (int i = n-1; i>0; i--){
        v[i] = (u[i] - v[i+1]*c[i])/d[i];
    }

};

void LUdecomposition::print_solution(){

    cout << "solution v = [";
    for (int i = 0; i<n; i++){
        cout << " " << v[i];
    }
    cout << " ]" << endl;
}