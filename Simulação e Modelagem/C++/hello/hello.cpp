#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

vec smooth_interpolate(double t, double tf, double x0, double xf){
	double dx = xf - x0;
	double t1 = t/tf;
	double t2 = t1*t1;
	double t3 = t2*t1;
	double t4 = t3*t1;
	double t5 = t4*t1;
    return {x0 + dx*(10*t3 - 15*t4 + 6*t5), (30*dx/tf)*(t2 - 2*t3 + t4), (60*dx/(tf*tf))*(t1 - 3*t2 + 2*t3)};
}

 //typedef vec (*func)(double t);
 //func smooth2(double tf, double x0, double xf){
 //    func f3 = [] (double t) -> vec { return smooth_interpolate(t,1,0,1); };
 //    return f3;
 //}

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
    cout << smooth_interpolate(0.25, 1, 0, 1) << endl;

    typedef int (*func2)(int a);
    func2 f2 = [] (int a) -> int { return a; };
    cout << f2(3) << endl;

    typedef vec (*func3)(double t);
    func3 f3 = [] (double t) -> vec { return smooth_interpolate(t,1,0,1); };
    cout << f3(0.25) << endl;

    return 0;
}