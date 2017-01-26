#include "Utilities.h"

#define SIM_TIME 4000

int main(){
	uint32_t t1;
	uint32_t t2;
	int deltaT=0;
	int T=1; 

	Setup();
	encSetup();
	
	long x_ [SIM_TIME];
	int i = 0;

	V2pwm(0.0);  

	t1=millis();
  
	deltaT=millis()-t1;
    while(deltaT<2000){
		deltaT=millis()-t1;
	}
	  t1=millis();
	  t2=millis();
	  V2pwm(6.0);
	  
  while(millis()-t2 < SIM_TIME){       
	  
	  LeEncoder();
		 
       deltaT=millis()-t1;
       if(deltaT>=T){
		   t1=millis();
           x_[i] = ReturnCounter();
           i++;
	   }
  }
  V2pwm(0.0);

  for (i=0; i<SIM_TIME; i++){
	  printf("%ld\n", x_[i]);
  }
return 0;

}
