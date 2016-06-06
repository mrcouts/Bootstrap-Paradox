#include <cmath>
#include <algorithm>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

class RK {
public:
    RK(string method);
    ~RK();
    void Doit();

    int N;
    vec *a__;
    vec b_;
    vec c_;
    
   };