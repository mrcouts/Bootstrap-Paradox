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

    vec *a;
    vec b;
    vec c;
    int N;
   };