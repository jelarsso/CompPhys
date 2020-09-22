#ifndef JACOBI_ROTATION_H
#define JACOBI_ROTATION_H


int jacobi_rotate(double** array, double* eigvals, int size, double tolerance);
double off_diagonl(double** array, int size);
void rotate(double** array, int size, int l, int k, double c, double s);

#endif