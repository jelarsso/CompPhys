#include "algorithms.hpp"
#include<armadillo>
#include<iostream>
#include<string>
#include<iomanip>

int inline periodic(int indx, int size){
    /*
    Simple periodic function to allow periodic boundaries.
    */
    return (indx+size)%size;
};

std::ofstream output_file;
void output(int n, arma::Col<double> *dump, std::string filename){
    /*
    output the one dimensional data in the vector dump to filename.
    */
    if (output_file.is_open()==false){
        output_file.open(filename);
    }
    for (int i=0;i<n+1;i++){
        output_file << (*dump)(i) << " ";
    }
    output_file << "\n";
};

void output(int n, arma::Mat<double> *dump, std::string filename){
    /*
    output the (n x n) two dimensional data in the matrix dump to filename.
    */
    if (output_file.is_open()==false){
        output_file.open(filename);
    }
    for (int i=0;i<n+1;i++){
        for (int j=0; j<n+1;j++) output_file << (*dump)(i,j) << " ";
        output_file << "\n";
    }
    output_file << "\n";
};

void output(int nx, int ny, arma::Mat<double> *dump, std::string filename){
    /*
    output the (nx x ny) two dimensional data in the matrix dump to filename.
    */

    if (output_file.is_open()==false){
        output_file.open(filename);
    }
    
    for (int i=0;i<nx+1;i++){
        for (int j=0; j<ny+1;j++) output_file << std::setprecision(15) << (*dump)(i,j) << " ";
        output_file << "\n";
    }
    output_file << "\n";
};


void close(){
    output_file.close(); // to make sure the file is closed
}

void lusolve(double a, double b, double c, int n, arma::Col<double>* rhs, arma::Col<double>* sol){
    //LUsolve from project 1, only used by the functions below.
    // solution is in sol;

    n = n-2;
    arma::Col<double> l(n),d(n+1), u(n+1);
    d(0) = b;
    for (int i = 1; i<n+1; i++){
        l(i-1) = a/d(i-1);
        d(i) = b - l(i-1)*c;
    }

    u(0) = (*rhs)(1);
    for (int i = 1; i<n+1; i++){
        u(i) = (*rhs)(i+1) - u(i-1)*l(i-1);
    }
    
    (*sol)(n+1) = u(n)/d(n);
    for (int i = n-1; i>=0; i--){
        (*sol)(i+1) = (u(i) - (*sol)(i+2)*c)/d(i);
    }
};

void forward_euler(int n, double dx, double alpha, int number_of_steps_t, arma::Col<double> *usolve,std::string filename){
    /*
    perfrom the forward euler in one dimension with int n + 1 grid points an dx grid spacings. alpha = dt/dx/dx and int number_of_steps_t is self-explanantory.
    *usolve is a pointer to a armadillo column vector of size n+1 where the solution is stored.
    The data is written to filename during the simulation.
    */
    arma::Col<double> u(n+1);
    for (int i=1; i<n; i++){
        u(i) = init_func(dx*i);
    }
    u(0) = 0;
    u(n) = 0;

    output(n,&u,filename);

    for (int t=0; t<number_of_steps_t;t++){
        for (int i=1; i<n; i++){
            (*usolve)(i) = alpha*u(i-1) + (1-2*alpha)*u(i) + alpha*u(i+1);
        }
        (*usolve)(0) = 0;
        (*usolve)(n) = 0;
        for (int j = 0; j<n+1;j++) u(j) = (*usolve)(j);
    output(n,usolve,filename);
    }
    output_file.close();
};

void backward_euler(int n, double dx, double alpha, int number_of_steps_t, arma::Col<double> *usolve,std::string filename){
    /*
    perfrom the backward euler in one dimension with int n + 1 grid points an dx grid spacings. alpha = dt/dx/dx and int number_of_steps_t is self-explanantory.
    *usolve is a pointer to a armadillo column vector of size n+1 where the solution is stored.
    The data is written to filename during the simulation.
    */
    arma::Col<double> uprev(n+1);
    for (int i=1; i<n; i++){
        (*usolve)(i) = init_func(dx*i);
        uprev(i) = init_func(dx*i);
    }
    (*usolve)(0) = 0;
    uprev(0) = 0;
    (*usolve)(n) = 0;
    uprev(n) = 0;

    output(n,usolve,filename);

    
    for (int t = 0; t < number_of_steps_t; t++){
        lusolve(-alpha,1+2*alpha,-alpha, n, &uprev, usolve);
        (*usolve)(0) = 0;
        (*usolve)(n) = 0;
        for (int j = 0; j<n+1;j++) uprev(j) = (*usolve)(j);
    output(n,usolve,filename);
    }
    output_file.close();
};


void cranky_nicholson(int n, double dx, double alpha, int number_of_steps_t,arma::Col<double> *usolve, std::string filename){
    /*
    perfrom the Crank Nicolson algorithm in one dimension with int n + 1 grid points an dx grid spacings. alpha = dt/dx/dx and int number_of_steps_t is self-explanantory.
    *usolve is a pointer to a armadillo column vector of size n+1 where the solution is stored.
    The data is written to filename during the simulation.
    */
    arma::Col<double> rhs(n+1);

    for (int i=1; i<n; i++){
        (*usolve)(i) = init_func(dx*i);
    }
    (*usolve)(0) = 0;
    (*usolve)(n) = 0;
    output(n,usolve,filename);

    for (int t = 0; t<number_of_steps_t; t++){
        for (int i = 1; i<n;i++) rhs(i) = alpha*(*usolve)(i-1) + (2-2*alpha)*(*usolve)(i)  + alpha*(*usolve)(i+1);
        rhs(0) = 0;
        rhs(n) = 0;
        lusolve(-alpha,2+2*alpha,-alpha,n,&rhs,usolve);
        (*usolve)(0)=0;
        (*usolve)(n)=0;
        output(n,usolve,filename);
    }
    output_file.close();
    
};

void forward_euler2d(int n, double dx, double alpha, int number_of_steps_t, arma::Mat<double> *usolve,std::string filename){
    /*
    perfrom the forward euler in two dimension with int (n + 1)x(n +1) grid points an dx grid spacings. alpha = dt/dx/dx and int number_of_steps_t is self-explanantory.
    *usolve is a pointer to a armadillo matrix of size (n+1)x(n+1) where the solution is stored.
    The data is written to filename during the simulation.
    */
    // boundaries are 0.
    arma::Mat<double> u(n+1,n+1,arma::fill::zeros);
    for (int i=1; i<n; i++){
        for (int j=1;j<n;j++){
            u(i,j) = init_func(dx*i,dx*j);
            (*usolve)(i,j) = init_func(dx*i,dx*j);
            }
    }
    
    output(n,usolve,filename);
    
    for (int t=0; t<number_of_steps_t;t++){
        for (int i=1; i<n; i++){
            for (int j=1;j<n;j++) (*usolve)(i,j) = u(i,j) + alpha*(u(i-1,j) + u(i+1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j));
        }
        for (int ix = 0; ix<n+1;ix++){
            for (int iy=0;iy<n+1;iy++) u(ix,iy) = (*usolve)(ix,iy);
        }
    output(n,usolve,filename);
    }
    output_file.close();
};

// The following two functions are specific adaptions of forward_euler2d to the lithosphere simulations and can probably not be understood without that context provided by the paper and project description.

void forward_euler2d_litho(int nx, int ny, double dx, double alpha, double boundary_upper, double boundary_lower, double qupper, double qmiddle,double qlower, double w_subduct,int number_of_steps_t, arma::Mat<double> *usolve, std::string filename){
    /* nb: all variables must be scaled beforehand
    perfrom the forward euler in two dimension with int (nx + 1)x(ny +1) grid points an dx grid spacings.
    This allows specific boundary conditions in along the y-edges (x=0 or x=nx*dx) set by the inital condition read from *usolve.
    qupper,qlower,qmiddle is the heat production in the three zones described in the method section of the paper, the division between these three region is set by boundary_lower/upper.
    w_subduct is the width of the subduction area that is enriched.
    alpha = dt/dx/dx and int number_of_steps_t is self-explanantory.
    *usolve is a pointer to a armadillo matrix of size (n+1)x(n+1) where the intial condition is stored beforehand, after the simulation the solution is stored here.
    The data is written to filename during the simulation.
    */
    arma::Mat<double> u(nx+1,ny+1,arma::fill::zeros);

    int i_bl = (int) (boundary_lower/dx);
    int i_bu = (int) (boundary_upper/dx);
    double dt = alpha*dx*dx;

    int half_subduct = (int) (w_subduct/dx/2);
    
    for (int i=0; i<nx+1; i++){
        for (int j=0;j<ny+1;j++){
            u(i,j) = (*usolve)(i,j);
        }
    }
    
    
    output(nx,ny,usolve,filename);
    
    for (int t=0; t<number_of_steps_t;t++){
        for (int i=1; i<nx; i++){
            for (int j=1;j<ny;j++){
                (*usolve)(i,j) = u(i,j) + alpha*(u(i-1,j) + u(i+1,j) + u(i,j+1) + u(i,j-1) - 4*u(i,j));
                if (i<i_bu){
                    (*usolve)(i,j) += qupper*dt;
                }else if (i<i_bl){
                    (*usolve)(i,j) += qmiddle*dt;
                }else{
                    if (j> ny/2 - half_subduct && j< ny/2 +half_subduct){
                        (*usolve)(i,j) += qlower*dt + radiohalflife(t*dt)*dt;
                    }else{
                        (*usolve)(i,j) += qlower*dt;
                    }
                    
                }
            } 
        }
        for (int ix = 0; ix<nx+1;ix++){
            for (int iy=0;iy<ny+1;iy++) u(ix,iy) = (*usolve)(ix,iy);
        }
    if (t%100==0) output(nx,ny,usolve,filename);
    }
    output_file.close();
};


void forward_euler2d_litho_pb(int nx, int ny, double dx, double alpha, double boundary_upper, double boundary_lower, double qupper, double qmiddle,double qlower,int number_of_steps_t, arma::Mat<double> *usolve, std::string filename){
    /* nb: all variables must be scaled beforehand
    perfrom the forward euler in two dimension with int (nx + 1)x(ny +1) grid points an dx grid spacings.
    This allows specific boundary conditions in along the y-edges (x=0 or x=nx*dx) set by the inital condition read from *usolve and uses PERIODIC BOUNDARIES along the x-edges!!
    qupper,qlower,qmiddle is the heat production in the three zones described in the method section of the paper, the division between these three region is set by boundary_lower/upper.
    w_subduct is the width of the subduction area that is enriched.
    alpha = dt/dx/dx and int number_of_steps_t is self-explanantory.
    *usolve is a pointer to a armadillo matrix of size (n+1)x(n+1) where the intial condition is stored beforehand, after the simulation the solution is stored here.
    The data is written to filename during the simulation.
    */
    arma::Mat<double> u(nx+1,ny+1,arma::fill::zeros);
    double dt = alpha*dx*dx;

    int i_bl = (int) (boundary_lower/dx);
    int i_bu = (int) (boundary_upper/dx);

    for (int i=0; i<nx+1; i++){
        for (int j=0;j<ny+1;j++){
            u(i,j) = (*usolve)(i,j);
        }
    }
    
    output(nx,ny,usolve,filename);
    
    for (int t=0; t<number_of_steps_t;t++){
        for (int i=1; i<nx; i++){
            for (int j=0;j<ny+1;j++){
                (*usolve)(i,j) = u(i,j) + alpha*(u(i,periodic(j-1,ny)) + u(i,periodic(j+1,ny)) + u(i+1,j) + u(i-1,j) - 4*u(i,j));
                if (i<i_bu){
                    (*usolve)(i,j) += qupper*dt;
                }else if (i<i_bl){
                    (*usolve)(i,j) += qmiddle*dt;
                }else{
                    (*usolve)(i,j) += qlower*dt;
                }
            } 
        }
        for (int ix = 0; ix<nx+1;ix++){
            for (int iy=0;iy<ny+1;iy++) u(ix,iy) = (*usolve)(ix,iy);
        }
    if (t%100==0) output(nx,ny,usolve,filename);
    }
    output_file.close();
};