#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

mat join_diag(mat A, mat B){
    return join_vert( join_horiz(A, zeros(A.n_rows, B.n_cols)), join_horiz(zeros(B.n_rows, A.n_cols), B) );
}

int main(void){
	cout << "Hello! This is a C++ program." << endl;
	vec A = randn(4);
    vec B = randn(2);
    
    field<vec> F(2);
    F(0) = A;
    F(1) = B; 
    
    //F.print("F:");
    cout << F << endl;
    cout << F(0) << endl;
    cout << F(1) << endl;

    mat C = randn(4,4);
    mat D = randn(2,2);
    cout << C << endl;
    cout << D << endl;
    cout << join_diag( join_diag(C,D), C ) << endl;

    return 0; }