#include <cmath>
#include <algorithm>
#include <armadillo>
#include "SomeUtilities.h"

using namespace std;
using namespace arma;

class Serial {
public:
    Serial(int dof, vec l_, vec lg_, vec m_, cube I__, mat (*fDH)(vec, vec, vec));
    ~Serial();
    void Doit(vec q0_);

    int dof;
    vec l_;
    vec lg_;
    vec m_;
    cube I__;
    cube Hr__;
    cube H__;
    cube z__;
    cube o__;
    cube og__;
    cube Jv__;
    cube Jw__;
    cube Jw2__;
    mat Jv_n_;
    mat Jw_n_;
    mat M_;
    mat (*fDH)(vec, vec, vec); };