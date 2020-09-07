#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <string>


using namespace std;
using namespace arma;

int main(int argv, char* argc[]){

  int n = atoi(argc[1]);
  double h = 1/double(n-1);
  mat A = Mat<double>(n,n,fill::zeros);
  vec b = vector<double>(n);


  for (int i=0;i<n;i++){
    A(i,i) = 2;
  }

  for (int i=0;i<n-1;i++){
    A(i,i+1) = -1;
    A(i+1,i) = -1;
  }


  for (int i=0;i<n;i++){
      double x = double(i)*h;
      b(i) = 100*exp(-10*x)*h*h;
  }


  mat L,U;

  lu(L,U,A);

  /*
  for (int i = 0; i<n; i++){
    for (int j = 0; j<n; j++){

      cout << A(i,j) << " ";
    }
    cout << endl;
  }*/


  vec u = solve(L,b);
  
  vec v = solve(U,u);

  //ofstream file_to_write;
  //string filename = "armadillo_output.data";
  //file_to_write.open(filename);

  //for (int i=0;i<n;i++){
  //  file_to_write << v(i) << "\n";
  //}
  //file_to_write.close();
  

  return 0;
}
