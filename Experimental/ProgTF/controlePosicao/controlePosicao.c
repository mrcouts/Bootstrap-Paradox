#include "Utilities.h"

int main(){
	uint32_t t1;
	uint32_t t2;
	int deltaT=0;
	int deltaT2=0;
	int T=1;
	int T2=2000;
	
	Setup();
	encSetup();
	
	float ref = 1;// 2*PI;
	float ref2 = 0.0;
	float ref2_1 = 0.0;
	float theta = 0.0;
	float e = 0.0;
	float e_1 = 0.0;
	float v = 0.0;
	float v_1 = 0.0;
	
	float a1 = 2105.89;
	float a2 = -2064.47;
	float a3 = 0.399521;

	V2pwm(0.0);  
	t1=millis();
	t2=millis();
	  
  while(1){
	   LeEncoder();		 
       deltaT=millis()-t1;
       deltaT2=millis()-t2;
       if(deltaT>=T){
		   t1=millis();
		   theta = -0.0001*PI*ReturnCounter();
		   //ref2 = 0.99203187251*ref2_1 + 0.00796812749*ref;
		   ref2 = 0.9950124688*ref2_1 + 0.004987531*ref;	
		   e = ref2 - theta;
		   v = -a3*v_1 + a1*e + a2*e_1;
          V2pwm( v );
          v_1 = v;
          e_1 = e;
          ref2_1 = ref2;  
	   }
	   
	   if(deltaT2>=T2){
		   t2=millis();
		   ref = -ref;
	   }
  }

return 0;
}
