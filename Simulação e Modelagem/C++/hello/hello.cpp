#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(void){
	cout << "Hello! This is a C++ program." << endl;
	vec A = randn(4);
    vec B = randn(2);
    vec B2 = randn(3);
    
    field<vec> F(3);
    F(0) = A;
    F(1) = B;
    F(2) = B2;

    //F.print("F:");
    cout << F << endl;
    cout << F(0) << endl;
    cout << F(1) << endl;

    mat C;
    C << 1 << 2 << 3 << 4 << 5 << 6 << endr
      << 0 << 0 << 0 << 0 << 0 << 0 << endr
      << 0 << 0 << 0 << 0 << 0 << 0 << endr
      << 0 << 0 << 0 << 0 << 0 << 0 << endr;
    cout << C << endl;
    cout << arma::rank(C) << endl;
    cout << C.col(0) << endl;
    cout << C.cols(2,5) << endl;
    mat C2 = join_horiz(C.col(0),C.cols(2,5));
    cout << C2 << endl;
    return 0;
}