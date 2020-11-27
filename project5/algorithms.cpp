#include "algorithms.hpp"
#include<armadillo>
#include<iostream>
#include<string>


std::ofstream output_file;
void output(int n, arma::Col<double> *dump, std::string filename){
    if (output_file.is_open()==false){
        output_file.open(filename);
    }
    for (int i=0;i<n+1;i++){
        output_file << (*dump)(i) << " ";
    }
    output_file << "\n";
};

void output(int n, arma::Mat<double> *dump, std::string filename){
    if (output_file.is_open()==false){
        output_file.open(filename);
    }
    for (int i=0;i<n+1;i++){
        for (int j=0; j<n+1;j++) output_file << (*dump)(i) << " ";
        output_file << "\n";
    }
    output_file << "\n";
};



void lusolve(double a, double b, double c, int n, arma::Col<double>* rhs, arma::Col<double>* sol){
    // solution is in sol
    arma::Col<double> l(n),d(n+1), u(n+1);
    d(0) = b;
    for (int i = 1; i<n+1; i++){
        l(i-1) = a/d(i-1);
        d(i) = b - l(i-1)*c;
    }

    u(0) = (*rhs)(0);
    for (int i = 1; i<n; i++){
        u(i) = (*rhs)(i) - u(i-1)*l(i-1);
    }

    
    u(n) = u(n)/d(n);
    (*sol)(n-1)  = u(n);
    for (int i = n-2; i>=1; i--){
        (*sol)(i) = (u(i) - (*sol)(i+1)*c)/d(i);
    }
};

void forward_euler(int n, double dx, double alpha, int number_of_steps_t, arma::Col<double> *usolve,std::string filename){
    arma::Col<double> u(n+1);
    for (int i=1; i<n; i++){
        u(i) = init_func(dx*i);
    }
    u(0) = 0;
    u(n) = 0;
    
    for (int t=0; t<number_of_steps_t;t++){
        for (int i=1; i<n; i++){
            (*usolve)(i) = alpha*u(i-1) + (1-2*alpha)*u(i) + alpha*u(i+1);
        }
        for (int j = 0; j<n+1;j++) u(j) = (*usolve)(j);
    output(n,usolve,filename);
    }
    output_file.close();
};

void backward_euler(int n, double dx,double alpha, int number_of_steps_t, arma::Col<double> *usolve,std::string filename){
    arma::Col<double> uprev(n+1);
    (*usolve)(0) = 0;
    uprev(0) = 0;
    (*usolve)(n) = 0;
    uprev(n) = 0;
    for (int i=1; i<n; i++){
        (*usolve)(i) = init_func(dx*i);
        uprev(i) = init_func(dx*i);
    }
    
    for (int t = 0; t < number_of_steps_t; t++){
        lusolve(-alpha,1+2*alpha,-alpha, n+1, &uprev, usolve);
        (*usolve)(0) = 0;
        (*usolve)(n) = 0;
        for (int j = 0; j<n+1;j++) uprev(j) = (*usolve)(j);
    output(n,usolve,filename);
    }
    output_file.close();
};


void cranky_nicholson(int n, double dx, double alpha, int number_of_steps_t,arma::Col<double> *usolve, std::string filename){
    arma::Col<double> rhs(n+1);
    (*usolve)(0) = 0;
    (*usolve)(n) = 0;

    for (int i=1; i<n; i++){
        (*usolve)(i) = init_func(dx*i);
    }

    for (int t = 0; t<number_of_steps_t; t++){
        for (int i = 1; i<n;i++) rhs(i) = alpha*(*usolve)(i-1) + (2-2*alpha)*(*usolve)(i)  + alpha*(*usolve)(i+1);
        rhs(0) = 0;
        rhs(n) = 0;
        lusolve(-alpha,2+2*alpha,-alpha,n,&rhs,usolve);
        (*usolve)(0)=0;
        (*usolve)(n)=0;
        output(n,usolve,filename);
    }
    
};

void forward_euler2d(int n, double dx, double alpha, int number_of_steps_t, arma::Mat<double> *usolve,std::string filename){
    arma::Mat<double> u(n+1,n+1);
    for (int i=1; i<n; i++){
        for (int j=1;j<n;j++) u(i) = init_func(dx*i,dx*j);
    }
    for (int j=0;j<n+1;j++){
    u(0,j) = 0;
    u(n,j) = 0;
    u(j,0) = 0;
    u(j,n) = 0;
    }
    
    for (int t=0; t<number_of_steps_t;t++){
        for (int i=1; i<n; i++){
            for (int j=1;j<n;j++) (*usolve)(i,j) = u(i,j) + alpha*(u(i-1,j) + u(i+1,j) + u(i,j+1) + u(i,j-1) + 4*u(i,j));
        }
        for (int ix = 0; ix<n+1;ix++){
            for (int iy=0;iy<n+1;iy++) u(ix,iy) = (*usolve)(ix,iy);
        }
    output(n,usolve,filename);
    }
    output_file.close();
};