#include <cmath>
#include <algorithm>
#include <armadillo>
#include "SomeUtilities.h"

using namespace std;
using namespace arma;

class Serial {
public:
    Serial(int dof, vec l_, vec lg_, mat (*fDH)(vec, vec, vec));
    ~Serial();
    void Doit(vec q0_);

    int dof;
    vec l_;
    vec lg_;
    cube Hr__;
    cube H__;
    mat (*fDH)(vec, vec, vec); };