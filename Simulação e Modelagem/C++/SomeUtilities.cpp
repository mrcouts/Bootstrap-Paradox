#include "SomeUtilities.h"

mat Rotx(double theta){
	mat M;
	M << 1 << 0          <<  0          << endr
      << 0 << cos(theta) << -sin(theta) << endr
      << 0 << sin(theta) <<  cos(theta) << endr;
    return M;}

mat Roty(double theta){
	mat M;
	M <<  cos(theta) << 0 << sin(theta) << endr
      <<  0          << 1 << 0          << endr
      << -sin(theta) << 0 << cos(theta) << endr;
    return M;}

mat Rotz(double theta){
	mat M;
	M << cos(theta) << -sin(theta) << 0 << endr
      << sin(theta) <<  cos(theta) << 0 << endr
      << 0          << 0           << 1 << endr;
    return M;}

mat Hx(double theta, double dx, double dy, double dz){
	mat M = Rotx(theta);
	mat v1;
	v1 << dx << endr
	   << dy << endr
	   << dz << endr;
	mat v2;
	v2 << 0 << 0 << 0 << 1 << endr;
	return join_vert(join_horiz(M,v1),v2);}

mat Hy(double theta, double dx, double dy, double dz){
	mat M = Roty(theta);
	mat v1;
	v1 << dx << endr
	   << dy << endr
	   << dz << endr;
	mat v2;
	v2 << 0 << 0 << 0 << 1 << endr;
	return join_vert(join_horiz(M,v1),v2);}

mat Hz(double theta, double dx, double dy, double dz){
	mat M = Rotz(theta);
	mat v1;
	v1 << dx << endr
	   << dy << endr
	   << dz << endr;
	mat v2;
	v2 << 0 << 0 << 0 << 1 << endr;
	return join_vert(join_horiz(M,v1),v2);}

mat Ht(double dx, double dy, double dz){
	mat M = eye<mat>(3,3);
	mat v1;
	v1 << dx << endr
	   << dy << endr
	   << dz << endr;
	mat v2;
	v2 << 0 << 0 << 0 << 1 << endr;
	return join_vert(join_horiz(M,v1),v2);}

mat H(char axis, double theta, double dx, double dy, double dz){
	if(axis == 'x')
		return Hx(theta, dx, dy, dz);
	else if(axis == 'y')
		return Hy(theta, dx, dy, dz);
	else if(axis == 'z')
		return Hz(theta, dx, dy, dz);
	else
		return Ht(dx, dy, dz);}

mat H_d(double a, double alpha, double d, double theta){
	return Hz(theta, 0, 0, 0)*Hx(alpha, a,0,d);}