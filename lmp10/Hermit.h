#ifndef _HERMIT_H_
#define _HERMIT_H

#include "piv_de_solver.h"
#include "makespl.h"

double H(int step , double x);

double gdy_0(double a , int i , int k);

double gdy_1(double a , int i , int k);

double gdy_2(double a , int i , int k);

double gdy_3(double a , int i , int k);

double dla_y(double a);


#endif
