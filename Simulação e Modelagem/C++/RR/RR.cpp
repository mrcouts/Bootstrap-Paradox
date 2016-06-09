#include "RR.h"

mat fDH_RR(vec q0_, vec l_, vec lg_){
    mat M;
    M << l_(0) << 0 << 0 << q0_(0) << -l_(0) + lg_(0) << 0 << 0 << true << endr
      << l_(1) << 0 << 0 << q0_(1) << -l_(1) + lg_(1) << 0 << 0 << true << endr;
    return M;}