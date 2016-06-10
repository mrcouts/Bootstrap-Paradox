#include "Dy.h"

Dy::Dy(int dof){
	Mh_.eye(dof,dof);
	vh_.zeros(dof);
	gh_.zeros(dof); }

Dy::~Dy(){
	Mh_.clear();
	vh_.clear();
	gh_.clear(); }