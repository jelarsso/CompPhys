#include <iostream>
#include <vector>
#include <armadillo>


using namespace std;
using namespace arma;

int main(int argv, char* argc[]){


  int n = atoi(argc[0]);
  mat A = Mat<double>(n,n,fill::zeros);
  vec b = vector<double>(n);

  mat L,U;

  lu(L,U,A);

  vec u = solve(L,b);
  vec v = solve(U,u);


  return 0;
}
