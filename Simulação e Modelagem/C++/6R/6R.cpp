#include "6R.h"

mat fDH_6R(vec q0_, vec l_, vec lg_){
    mat M;
    M << 0 << PI/2 << 0             << q0_(0) + PI/2 << 0 << 0               << lg_(0)          << true << endr
      << 0 << PI/2 << l_(0) + l_(1) << q0_(1) + PI   << 0 << -l_(1) + lg_(1) << 0               << true << endr
      << 0 << PI/2 << 0             << q0_(2) + PI   << 0 << 0               << lg_(2)          << true << endr
      << 0 << PI/2 << l_(2) + l_(3) << q0_(3) + PI   << 0 << -l_(3) + lg_(3) << 0               << true << endr
      << 0 << PI/2 << 0             << q0_(4) + PI   << 0 << 0               << lg_(4)          << true << endr
      << 0 << 0    << l_(4) + l_(5) << q0_(5)        << 0 << 0               << -l_(5) + lg_(5) << true << endr;
    return M;}