#include "Filter.h"

Filter::Filter(int size, mat AB_){
    this->order = AB_.n_cols-1;
    this->size = size;
    this->a_ = (AB_.row(0)/AB_(0,0)).t();
    this->b_ = (AB_.row(1)/AB_(0,0)).t();
    u__.zeros(size,order+1);
    y__.zeros(size,order+1);
}

Filter::~Filter(){
    a_.clear();
    b_.clear();
    u__.clear();
    y__.clear();
}

vec Filter::Doit(vec u_){
    u__.col(0) = u_;
    y__.col(0) = b_(0)*u__.col(0);
    for(int i = 1; i <= order; i++)
        y__.col(0)  += b_(i)*u__.col(i) - a_(i)*y__.col(i);

    
    for(int i = order; i > 0; i--){
        u__.col(i) = u__.col(i-1);
        y__.col(i) = y__.col(i-1);
    }

    //u__ = shitf(u__,+1);
    //y__ = shitf(u__,+1);

    return y__.col(0);
}

mat Tustin(double T, double w0, mat ABs_){
    int order = ABs_.n_cols-1;
    mat AB_ = ABs_;
    double c;
    if (w0 == 0.0) c = 2/T;
    else c = w0/tan(0.5*w0*T);
    double c2 = c*c;
    double c3 = c2*c;
    double c4 = c3*c;
    double c5 = c4*c;

    switch(order){
        case 1:
            AB_.col(0) = ABs_.col(0) + c*ABs_.col(1);
            AB_.col(1) = ABs_.col(0) - c*ABs_.col(1);
            break;
        case 2:
            AB_.col(0) =   ABs_.col(0) + c*ABs_.col(1) + c2*ABs_.col(2);
            AB_.col(1) = 2*ABs_.col(0)               - 2*c2*ABs_.col(2);
            AB_.col(2) =   ABs_.col(0) - c*ABs_.col(1) + c2*ABs_.col(2);
            break;
        case 3:
            AB_.col(0) =   ABs_.col(0) + c*ABs_.col(1) + c2*ABs_.col(2) +   c3*ABs_.col(3);
            AB_.col(1) = 3*ABs_.col(0) + c*ABs_.col(1) - c2*ABs_.col(2) - 3*c3*ABs_.col(3);
            AB_.col(2) = 3*ABs_.col(0) - c*ABs_.col(1) - c2*ABs_.col(2) + 3*c3*ABs_.col(3);
            AB_.col(3) =   ABs_.col(0) - c*ABs_.col(1) + c2*ABs_.col(2) -   c3*ABs_.col(3);
            break;
        case 4:
            AB_.col(0) =   ABs_.col(0) +   c*ABs_.col(1)  +   c2*ABs_.col(2) + c3*ABs_.col(3) +   c4*ABs_.col(4);
            AB_.col(1) = 4*ABs_.col(0) + 2*c*ABs_.col(1)                   - 2*c3*ABs_.col(3) - 4*c4*ABs_.col(4);
            AB_.col(2) = 6*ABs_.col(0)                    - 2*c2*ABs_.col(2)                  + 6*c4*ABs_.col(4);
            AB_.col(3) = 4*ABs_.col(0) - 2*c*ABs_.col(1)                   + 2*c3*ABs_.col(3) - 4*c4*ABs_.col(4);
            AB_.col(4) =   ABs_.col(0) -   c*ABs_.col(1)  +   c2*ABs_.col(2) - c3*ABs_.col(3) +   c4*ABs_.col(4);
            break;
        case 5:
            AB_.col(0) =    ABs_.col(0) +   c*ABs_.col(1)  +   c2*ABs_.col(2) +   c3*ABs_.col(3) +   c4*ABs_.col(4) +    c5*ABs_.col(5);
            AB_.col(1) =  5*ABs_.col(0) + 3*c*ABs_.col(1)  +   c2*ABs_.col(2) -   c3*ABs_.col(3) - 3*c4*ABs_.col(4) -  5*c5*ABs_.col(5);
            AB_.col(2) = 10*ABs_.col(0) + 2*c*ABs_.col(1)  - 2*c2*ABs_.col(2) - 2*c3*ABs_.col(3) + 2*c4*ABs_.col(4) + 10*c5*ABs_.col(5);
            AB_.col(3) = 10*ABs_.col(0) - 2*c*ABs_.col(1)  - 2*c2*ABs_.col(2) + 2*c3*ABs_.col(3) + 2*c4*ABs_.col(4) - 10*c5*ABs_.col(5);
            AB_.col(4) =  5*ABs_.col(0) - 3*c*ABs_.col(1)  +   c2*ABs_.col(2) +   c3*ABs_.col(3) - 3*c4*ABs_.col(4) +  5*c5*ABs_.col(5);
            AB_.col(5) =    ABs_.col(0) -   c*ABs_.col(1)  +   c2*ABs_.col(2) -   c3*ABs_.col(3) +   c4*ABs_.col(4) -    c5*ABs_.col(5);
            break;
    }
    AB_ = AB_/AB_(0,0);
    return AB_;
}

mat Bessel(double w, int order){
    mat ABs_ = ones(2, order+1);
    switch(order){
        case 1:
            ABs_(0,0) = w;
            ABs_(0,1) = 1.0;
            ABs_(1,0) = ABs_(0,0);
            ABs_(1,1) = 0.0;
            break;
        case 2:
            ABs_(0,0) = 1.61803*pow(w,2);
            ABs_(0,1) = 2.20320*w;
            ABs_(0,2) = 1.0;
            ABs_(1,0) = ABs_(0,0);
            ABs_(1,1) = 0.0;
            ABs_(1,2) = 0.0;
            break;
        case 3:
            ABs_(0,0) = 2.77179*pow(w,3);
            ABs_(0,1) = 4.86636*pow(w,2);
            ABs_(0,2) = 3.41749*w;
            ABs_(0,3) = 1.0;
            ABs_(1,0) = ABs_(0,0);
            ABs_(1,1) = 0.0;
            ABs_(1,2) = 0.0;
            ABs_(1,3) = 0.0;
            break;
        case 4:
            ABs_(0,0) = 5.25820*pow(w,4);
            ABs_(0,1) = 11.1154*pow(w,3);
            ABs_(0,2) = 10.0702*pow(w,2);
            ABs_(0,3) = 4.73055*w;
            ABs_(0,4) = 1.0;
            ABs_(1,0) = ABs_(0,0);
            ABs_(1,1) = 0.0;
            ABs_(1,2) = 0.0;
            ABs_(1,3) = 0.0;
            ABs_(1,4) = 0.0;
            break;
        case 5:
            ABs_(0,0) = 11.2128*pow(w,5);
            ABs_(0,1) = 27.2182*pow(w,4);
            ABs_(0,2) = 29.3643*pow(w,3);
            ABs_(0,3) = 17.8198*pow(w,2);
            ABs_(0,4) = 6.17942*w;
            ABs_(0,5) = 1.0;
            ABs_(1,0) = ABs_(0,0);
            ABs_(1,1) = 0.0;
            ABs_(1,2) = 0.0;
            ABs_(1,3) = 0.0;
            ABs_(1,4) = 0.0;
            ABs_(1,5) = 0.0;
            break;
    }
    return ABs_;
}

mat dBessel(double w, int order){
    mat AB_ = Bessel(w,order);
    AB_(1,1) = AB_(0,0);
    AB_(1,0) = 0;
    return AB_;
}

mat Bessel_d(double T, double w, int order){
    return Tustin(T,w,Bessel(w,order));
}

mat dBessel_d(double T, double w, int order){
    return Tustin(T,w,dBessel(w,order));
}