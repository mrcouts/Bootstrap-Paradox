#include "RK.h"

RK::RK(string method){
	if(method == "Heun"){
		N = 2;
		a__ = new vec[N-1];
		a__[0] = {1.0};
		b_ = {0.5, 0.5}; }
	else if(method == "RK3"){
		N = 3;
		a__ = new vec[N-1];
		a__[0] = {0.5};
		a__[1] = {-1.0, 2.0};
		b_ = {1.0/6, 2.0/3, 1.0/6}; }
	else if(method == "RK4"){
		N = 4;
		a__ = new vec[N-1];
		a__[0] = {0.5};
		a__[1] = {0, 0.5};
		a__[2] = {0, 0, 1.0};
		b_ = {1.0/6, 1.0/3, 1.0/3, 1.0/6}; }
	else if(method == "RK6"){
		N = 7;
		a__ = new vec[N-1];
		a__[0] = {0.5};
		a__[1] = {2.0/9, 4.0/9};
		a__[2] = {7.0/36, 2.0/9, -1.0/12};
		a__[3] = {-35.0/144, -55.0/36, 35.0/48, 15.0/8};
		a__[4] = {-1.0/360, -11.0/36, -1.0/8, 0.5, 1.0/10};
		a__[5] = {-41.0/260, 22.0/13, 43.0/156, -118.0/39, 32.0/195, 80.0/39};
		b_ = {13.0/200, 0, 11.0/40, 11.0/40, 4.0/25, 4.0/25, 13.0/200}; }
	else
		cout << "Método indisponível!" << endl;

	c_.zeros(N-1);
	for(int i = 0; i<N-1; i++)
		c_(i) = sum(a__[i]); }

RK::~RK(){
	delete[] a__;
	b_.clear() ;
	c_.clear(); }

void RK::Doit(){}