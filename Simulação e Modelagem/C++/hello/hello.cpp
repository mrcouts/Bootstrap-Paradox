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

    double lg1 = 0.5*0.12;
    double lg2 = 0.5*0.15;
    double m1 = 0.143;
    double m2 = 0.171;
    double Jz1 = 171.6e-6;
    double Jz2 = 320.6e-6;
    double sigma = 0.3;

    uint nPI = 6;
    cube lg__; lg__.zeros(2, 1, (uint)pow(2,6));
    cube m__; m__.zeros(2, 1, (uint)pow(2,6));
    cube Jz__; Jz__.zeros(2, 1, (uint)pow(2,6));

    uint j = 0;
    for(int i1=-1; i1<2; i1 = i1 + 2){
        for(int i2=-1; i2<2; i2 = i2 + 2){
            for(int i3=-1; i3<2; i3 = i3 + 2){
                for(int i4=-1; i4<2; i4 = i4 + 2){
                    for(int i5=-1; i5<2; i5 = i5 + 2){
                        for(int i6=-1; i6<2; i6 = i6 + 2){
                            lg__(0,0,j) = lg1*(1+sigma*i1);
                            lg__(1,0,j) = lg2*(1+sigma*i2);
                            m__(0,0,j) =  m1*(1+sigma*i3);
                            m__(1,0,j) =  m2*(1+sigma*i4);
                            Jz__(0,0,j) = Jz1*(1+sigma*i5);
                            Jz__(1,0,j) = Jz2*(1+sigma*i6);
                            j++;
                        }
                    }
                }
            }
        }
    }

    cout << Jz__ << endl;

    vec v = {1,2};
    mat M = ones(2,2);
    mat Const = v.t()*M*v;
    cout << Const(0,0) << endl;
    cout << v.t()*v << endl;

    cout << sign(randn(10)) << endl;
    cout << randu(10) % sign(randn(10)) << endl;
    return 0;
}