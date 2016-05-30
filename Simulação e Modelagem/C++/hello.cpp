#include <iostream>
#include "SomeUtilities.h"

//Simple C++ program
int main(void)
{
    cout << "Hello! This is a C++ program. \n";

    vec v;
    v.zeros(3);
    v(0) = 1;
    cout << v << endl;

    mat M;
    M << 1 << 2 << endr
      << 3 << 4 << endr; 
    cout << M << endl;

    mat M2 = Rotx(0.1);

    cout << M2        << endl;
    cout << Roty(0.1) << endl;
    cout << Rotz(0.1) << endl;
    cout << H('x',0.1,1,2,3) << endl;
    cout << H('y',0.1,1,2,3) << endl;
    cout << H('z',0.1,1,2,3) << endl;
    cout << H('0',0.1,1,2,3) << endl;
    cout << H_d(1, 0.1, 2, 0.2) << endl;
    return 0;
}