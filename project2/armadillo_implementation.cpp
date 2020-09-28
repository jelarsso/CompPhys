#include<armadillo>
#include"jacobi_rotation.hpp"
#include<iostream>
#include<stdlib.h>
#include<fstream>


int main(int argc, char *argv[]){
    int size = std::atoi(argv[1]);
    arma::mat A(size,size);
    A.zeros();

    double h = 1/(1.0*size);

    for (int i=0;i<size;i++){
        A(i,i) = 2/(h*h);// + (i+1)*h*(i+1)*h;
        if (i<size-1){
            A(i+1,i) = -1/(h*h);
            A(i,i+1) = -1/(h*h);
        }
    }


    arma::vec eigvalues;
    arma::mat eigvectors;

    arma::eig_sym(eigvalues, eigvectors, A);

    std::ofstream outfile;
    outfile.open("arma_data.out");
    for (int i=0;i<size+1;i++){
        for(int j=0;j<size;j++){
            if (i==0){
                outfile << eigvalues(j) << " ";
            }else{
                outfile << eigvectors(i-1,j) << " ";
            }
        }
        outfile << "\n";
    }

    return 0;
}