#ifndef _HERMIT_H_
#define _HERMIT_H


double H(int step , double x);

void sprawdz(spline_t *a , matrix_t *b);

double gdy_0(double a , int i , int k , double *x);

double gdy_1(double a , int i , int k , double *x);

double gdy_2(double a , int i , int k , double *x);

double gdy_3(double a , int i , int k , double*x);

double dla_y( int i , int ILE , double *x , double *y);


#endif
