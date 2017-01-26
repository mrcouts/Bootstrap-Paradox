#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <wiringPi.h>
#include <stdbool.h>

#define PWM_RANGE 1000
#define VCC 23.4
#define K (PWM_RANGE/VCC)
#define PI 3.141592653589793

void NPpwm(int x);
void V2pwm(float V);
void Setup();
int encSetup();
void LeEncoder();
long ReturnCounter();
void timerIsr();



#endif
