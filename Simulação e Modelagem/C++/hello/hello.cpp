#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(void){
	cout << "Hello! This is a C++ program." << endl;
    mat A; A.zeros(2,2);
    A << 1 << 2 << endr
      << 3 << 4 << endr;
    vec b = {1, 2};
    vec x = solve(A,b);
    A.print("A = ");
    b.print("b = ");
    x.print("x = ");
    return 0;
}