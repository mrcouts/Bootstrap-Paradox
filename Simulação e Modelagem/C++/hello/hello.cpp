#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

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
    return 0; }