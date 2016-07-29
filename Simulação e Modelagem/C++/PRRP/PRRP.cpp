#include "PRRP.h"

mat fDH_PRRP(vec q0_, vec l_, vec lg_){
    mat M;
    M << 0     << PI   << q0_(0)         << PI/2          << 0               << 0 << 0                   << false << endr
      << l_(1) << 0    << 0              << q0_(1)        << -l_(1) + lg_(1) << 0 << 0                   << true  << endr
      << 0     << PI/2 << 0              << q0_(2) + PI/2 << 0               << 0 << lg_(2)              << true  << endr
      << 0     << 0    << l_(2) + q0_(3) << 0             << 0               << 0 << q0_(3)-l_(3)+lg_(3) << false << endr;	
    return M;}