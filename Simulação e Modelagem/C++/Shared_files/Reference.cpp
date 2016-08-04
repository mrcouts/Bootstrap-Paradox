#include "Reference.h"

Reference::Reference(double tf, vec x0_, vec xf_){
    this->tf = tf;
    this->x0_ = x0_;
    dx_ = xf_ - x0_;
    dof = x0_.n_rows;
    r_ = x0_;
    dr_.zeros(dof);
    d2r_.zeros(dof);
}
Reference::~Reference(){
    r_.clear();
    dr_.clear();
    d2r_.clear();
    x0_.clear();
    dx_.clear();
}
void Reference::Doit(double t){
    this->t = t;
    double t1 = t/tf;
    double t2 = t1*t1;
    double t3 = t2*t1;
    double t4 = t3*t1;
    double t5 = t4*t1;
    if(t < 0){
        r_ = x0_;
        dr_ = zeros(dof);
        d2r_ = zeros(dof);
    }
    else if(t > tf){
        r_ = x0_ + dx_;
        dr_ = zeros(dof);
        d2r_ = zeros(dof);
    }
    else{
        r_ = x0_ + dx_*(10*t3 - 15*t4 + 6*t5);
        dr_ = (30*dx_/tf)*(t2 - 2*t3 + t4);
        d2r_ = (60*dx_/(tf*tf))*(t1 - 3*t2 + 2*t3);
    }  
}