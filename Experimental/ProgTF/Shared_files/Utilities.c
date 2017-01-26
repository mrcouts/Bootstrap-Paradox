#include "Utilities.h"

int EN = 22; //pino 15, gpio 19 DIGITAL
int INA = 27; //pino 13, gpio 19 DIGITAL
int INB = 17; //pino 11, gpio 19 DIGITAL
int PWM = 18; //pino 12, gpio 18 PWM
int encA = 6; //pino 31, gpio 6 
int encB = 12; //pino 32, gpio 12

int EstadoAtual = 0;
int EstadoAnterior = 0;
long counter = 0; // contador para o estado a cada 1ms.
bool encAr;
bool encBr;
volatile bool flag = 0;


void NPpwm(int x){
	if(x>=0){
		digitalWrite(INA, HIGH);
		digitalWrite(INB, LOW);
		pwmWrite(PWM, x);
	}
	else{
		digitalWrite(INA, LOW);
		digitalWrite(INB, HIGH);
		pwmWrite(PWM, -x);
	}
	
}

void V2pwm(float V){
	if(V > VCC) V = VCC;
	if(V < -VCC) V = -VCC;
	NPpwm((int)(K*V));	
}

void Setup(){
	if(wiringPiSetupGpio()==-1)
    exit(1);
  
    pinMode(INA,OUTPUT);
    pinMode(INB,OUTPUT);
    pinMode(EN,OUTPUT);	
    
    pinMode(PWM,PWM_OUTPUT);
    pwmSetMode(PWM_MODE_MS);
    pwmSetClock(2);
	pwmSetRange(PWM_RANGE);
    
    digitalWrite(EN, HIGH);
}

int encSetup(){
	pinMode(encA, INPUT);
	pinMode(encB, INPUT);
	
	if(wiringPiISR(encB, INT_EDGE_BOTH, &timerIsr)<0 || wiringPiISR(encA, INT_EDGE_BOTH, &timerIsr)<0){
		fprintf(stderr,"unable to setup ISR: %s \n", strerror(errno));
		return 1;
	}
	
	return 0;
}

void timerIsr(){
	flag = 1;
}

void LeEncoder(){
	if (flag){
		  encAr = digitalRead(encA); // Le encoder A
		  encBr = digitalRead(encB); // Le encoder B

		  if (!encAr && !encBr) {
			EstadoAtual = 0;
		  }
		  else if (encAr && !encBr) {
			EstadoAtual = 1;
		  }
		  else if (!encAr && encBr) {
			EstadoAtual = 2;
		  }
		  else if (encAr && encBr) {
			EstadoAtual = 3;
		  }

		  // criação da memória do estadoAnterior e incremento da variavel counter
		  switch (EstadoAtual) {
			case 0:
			  if (EstadoAnterior == 1) counter--;
			  else if (EstadoAnterior == 2) counter++;
			  break;
			case 1:
			  if (EstadoAnterior == 3) counter--;
			  else if (EstadoAnterior == 0) counter++;
			  break;
			case 2:
			  if (EstadoAnterior == 0) counter--;
			  else if (EstadoAnterior == 3) counter++;
			  break;
			case 3:
			  if (EstadoAnterior == 2) counter--;
			  else if (EstadoAnterior == 1) counter++;
			  break;
		  }

		  EstadoAnterior = EstadoAtual;
		 }
}


long ReturnCounter(){
	return counter;
}
