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
	vec v1 = {dx,dy,dz};
	vec v2 = {0,0,0,1};
	return join_vert(join_horiz(M,v1),v2.t());}

mat Hy(double theta, double dx, double dy, double dz){
	mat M = Roty(theta);
	vec v1 = {dx,dy,dz};
	vec v2 = {0,0,0,1};
	return join_vert(join_horiz(M,v1),v2.t());}

mat Hz(double theta, double dx, double dy, double dz){
	mat M = Rotz(theta);
	vec v1 = {dx,dy,dz};
	vec v2 = {0,0,0,1};
	return join_vert(join_horiz(M,v1),v2.t());}

mat Ht(double dx, double dy, double dz){
	mat M = eye<mat>(3,3);
	vec v1 = {dx,dy,dz};
	vec v2 = {0,0,0,1};
	return join_vert(join_horiz(M,v1),v2.t());}

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

mat join_diag(mat A, mat B){
    return join_vert( join_horiz(A, zeros(A.n_rows, B.n_cols)), join_horiz(zeros(B.n_rows, A.n_cols), B) ); }

vec join_vert(vec **lista, uint n){
    vec aux = *lista[0];
    for(uint i = 1; i<n; i++)
        aux = join_vert(aux, *lista[i]);
    return aux;
}

mat join_vert(mat **lista, uint n){
    mat aux = *lista[0];
    for(uint i = 1; i<n; i++)
        aux = join_vert(aux, *lista[i]);
    return aux;
}

mat join_diag(mat **lista, uint n){
    mat aux = *lista[0];
    for(uint i = 1; i<n; i++)
        aux = join_diag(aux, *lista[i]);
    return aux;
}

vec join_vert(vec v, uint n){
    vec aux = v;
    for(uint i = 1; i<n; i++)
        aux = join_vert(aux, v);
    return aux;
}

mat join_vert(mat v, uint n){
    mat aux = v;
    for(uint i = 1; i<n; i++)
        aux = join_vert(aux, v);
    return aux;
}
