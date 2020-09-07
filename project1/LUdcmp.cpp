#include "LUdcmp.hpp"
#include<iostream>
#include<fstream>
#include<string>


using namespace std;



LUdecomposition::LUdecomposition(int N, double *f_array){
    /*
    int N: size of array
    double *f_array: pointer to the start of the array containg the values of f_i

        the construcor for the class where a = c = -1 and b = 2
    */
    n = N;
    h_stepsize = 1/double(n+1);

    constant_abc = true;
    solved = false;
    LU_found = false;
    
    a = new double;
    b = new double;
    c = new double;

    *a = -1;
    *b = 2;
    *c = -1;
    
    f = f_array;
};


LUdecomposition::LUdecomposition(int N, double *a_array, double *b_array, double *c_array, double *f_array){
    /*
    int N: size of array
    double *a,*b,*c: pointer to the start of the array containing the N-1,N,N-1 values of a,b,c 
    double *f_array: pointer to the start of the array containing the N values of f_i

        the construcor for the class where a,b,c can be varying.
    */
    n = N;
    h_stepsize = 1/double(n+1);

    constant_abc = false;


    solved = false;
    LU_found = false;
    
    a = a_array;

    b = b_array;

    c = c_array;

    f = f_array;
};


LUdecomposition::~LUdecomposition(){
    // class destructor. clears the dynamically allocated arrays.
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
    /*
        find the LU decomposition
    */
    if (constant_abc){
        findLU_constant_abc();
    }else{
        findLU_varying_abc();
    }
    LU_found = true;
};


void LUdecomposition::findLU_constant_abc(){ //not needed for solving the system
    l = new double[n-1];
    d = new double[n];

    d[0] = *b;

    for (int i = 1; i<n; i++){
        l[i-1] = *a/d[i-1];
        d[i] = *b - l[i-1]*(*c);
    }
};

void LUdecomposition::findLU_varying_abc(){
    /*
    solves the LU decomposition for the system and stores it in l and d.
    */
    l = new double[n-1];
    d = new double[n];

    d[0] = b[0];

    for (int i = 1; i<n; i++){
        l[i-1] = a[i-1]/d[i-1];
        d[i] = b[i] - l[i-1]*c[i-1];
    }
};


void LUdecomposition::print_LU(){
    /*
        print the LU decomposition
    */
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
    /*
    solve and store result in v 
    */
    
    if (constant_abc){
        solve_constant_abc();
    }else{
        solve_varying_abc();
    }
    solved = true;
};


void LUdecomposition::solve_constant_abc(){
    /*
        solve with the special algorithm
    */
    u = new double[n];
    v = new double[n];

    u[0] = f[0];

    for (int i = 1; i<n; i++){
        u[i] = f[i] + u[i-1]*(double(i)/double(i+1));
    }
    u[n-1] = u[n-1]/(double(n+1)/double(n));
    v[n-1] = u[n-1];
    for (int i=n-2; i>=1; i--){
        v[i] = (u[i] + v[i+1])/(double(i+2)/double(i+1));
    }
};

void LUdecomposition::solve_varying_abc(){
    /*
        solve with the general case
    */
    if (!LU_found){
        throw "LU-decomposition must be found before calling this function";
    }

    u = new double[n];
    v = new double[n];


    u[0] = f[0];

    for (int i = 1; i<n; i++){
        u[i] = f[i] - u[i-1]*l[i-1];
    }
    
    u[n-1] = u[n-1]/d[n-1];
    v[n-1]  = u[n-1];

    for (int i = n-2; i>=1; i--){
        v[i] = (u[i] - v[i+1]*c[i])/d[i];
    }
};

void LUdecomposition::print_solution(){
    /*
        print the solution to stdout
    */

    cout << "solution v = [";
    for (int i = 0; i<n; i++){
        cout << " " << v[i];
    }
    cout << " ]" << endl;
};

void LUdecomposition::write_v_to_file(string filename){
    /*
        write the solution to the file string filename
    */
    ofstream file_to_write;
    file_to_write.open(filename);
    
    for (int i = 0; i<n; i++){
        file_to_write << v[i] << "\n";
    }
    file_to_write.close();
};
