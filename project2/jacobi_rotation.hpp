#include<string>
#ifndef JACOBI_ROTATION_H
#define JACOBI_ROTATION_H

int jacobi_rotate(double** array, double* eigvals, double** eigvectors, int size, int maxiter, double tolerance);
double off_diagonl(double** array, int size);
void rotate(double** array, double** eigvectors, int size, int l, int k, double c, double s);
void write_to_file(std::string filename, double** eigvectors, double* eigvals, int iters, int size);
void analytical(double* analyt_eigval, double** analyt_eigvec, int N);

#endif