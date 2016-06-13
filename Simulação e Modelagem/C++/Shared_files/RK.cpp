#include "RK.h"
#include <iostream>

RK::RK(string method, vec (*f_)(double, vec)){
	caso = 1;
	this->f_ = f_;
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

RK::RK(string method, Acceleration *AC){
	caso = 2;
	this->AC = AC;
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
	c_.clear();

	t_.clear();
	y__.clear();
	u__.clear();
	k__.clear(); }

void RK::Doit(double h, double tf, vec y0_){
	int nt = int(tf/h);
	t_.zeros(nt+1);
	for(int i = 0; i<nt+1; i++)
		t_(i) = i*h;
	y__.zeros(y0_.n_rows, 1, nt+1);
	u__.zeros(y0_.n_rows/2, 1, nt+1);
	y__.slice(0) = y0_;
	k__.zeros(y0_.n_rows, 1, N);

	vec aux_; aux_ = zeros(y0_.n_rows);
	for(int i = 0; i<nt; i++){
		switch (caso){
		    case 1: k__.slice(0) = f_(t_(i),y__.slice(i)); break;
		    case 2:
		        field<vec> F2(2);
		        F2 = AC->f2_(t_(i),y__.slice(i));
		        u__.slice(i) = F2(1); 
		        k__.slice(0) = F2(0);
		        break; }
		for(int j = 1; j<N; j++){
			aux_.zeros();
			for(int k = 0; k<j; k++)
				aux_ += a__[j-1](k)*k__.slice(k);
			switch (caso){
		        case 1: k__.slice(j) = f_(t_(i) + c_(j-1)*h, y__.slice(i) + h*aux_); break;
		        case 2: k__.slice(j) = AC->f_(t_(i) + c_(j-1)*h, y__.slice(i) + h*aux_); break; }
		aux_.zeros();
		for(int i = 0; i<N; i++)
			aux_ += b_(i)*k__.slice(i);
		y__.slice(i+1) = y__.slice(i) + h*aux_; }}
	if(caso == 2)
		u__.slice(nt) = AC->f2_(t_(nt),y__.slice(nt))(1); }