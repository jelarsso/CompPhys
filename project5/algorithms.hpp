#include<armadillo>
#include<fstream>

#ifndef algorithms
#define algorithms

double init_func(double);
double init_func(double,double);
void backward_euler(int, double, double, int, arma::Col<double>*,std::string);
void forward_euler(int, double, double, int, arma::Col<double>*,std::string);
void forward_euler2d(int,double,double,int,arma::Mat<double>*,std::string);
void cranky_nicholson(int n, double dx, double alpha, int number_of_steps_t, arma::Col<double>*,std::string);
void lusolve(double, double, double, int, arma::Col<double>*, arma::Col<double>*);
void output(int, arma::Col<double>*, std::string);
void output(int, arma::Mat<double>*, std::string);
#endif