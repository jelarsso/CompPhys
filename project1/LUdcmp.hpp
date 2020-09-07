#ifndef LUDCMP_H
#define LUDCMP_H
#include<string>

class LUdecomposition{
private:
int n; //size of matrix
double h_stepsize;
double *a;
double *b;
double *c;
double *d;
double *l;
double *f;
double *u;
double *v;
void solve_constant_abc();
void solve_varying_abc();
void findLU_constant_abc();
void findLU_varying_abc();
bool LU_found;
bool solved;
bool constant_abc;

public:
//constructors for the two use cases general and specialized algorithm
LUdecomposition(int N, double *a_array, double *b_array, double *c_array, double *f_array);
LUdecomposition(int N, double *f_array);

//destructor
~LUdecomposition();

void findLU();

double *geta();

void solve();

void print_solution();
void write_v_to_file(std::string filename);
void print_LU();
};
#endif