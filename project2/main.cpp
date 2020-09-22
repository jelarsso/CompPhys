#include "jacobi_rotation.hpp"
#include<iostream>




int main(){
    int size = 5;
    double** array;
    array = new double*[size];
    for (int i=0;i<size;i++){
        array[i] = new double[size];
    }

    array[0][0] = 1;
    array[1][1] = 5;
    array[2][2] = 7;
    array[3][3] = 3;
    array[4][4] = 5;

    array[0][1] = 6;
    array[1][0] = 6;
    array[2][0] = 7;
    array[0][2] = 7;
    array[0][3] = 1;
    array[3][0] = 1;
    array[0][4] = 6;
    array[4][0] = 6;

    array[1][2] = 3;
    array[2][1] = 3;
    array[3][1] = 4;
    array[1][3] = 4;
    array[1][4] = 5;
    array[4][1] = 5;

    array[2][3] = 9;
    array[3][2] = 9;
    array[2][4] = 4;
    array[4][2] = 4;

    array[4][3] = 3;
    array[3][4] = 3;

    for (int i = 0;i<size;i++){
        for (int j = 0;j<size;j++){
            std::cout << array[i][j] << " ";
        }
        std::cout << "\n";
    }

    int iters;
    double* eigvals;
    eigvals = new double[size];

    iters = jacobi_rotate(array,eigvals,size,1e-13);

    std::cout << "iters: " << iters << std::endl;
    std::cout << "eigenvalues" << std::endl;
    for (int i = 0; i<size; i++){
        std::cout << eigvals[i] << std::endl;
    }
    return 0;
};