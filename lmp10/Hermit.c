#include <stdio.h>
#include <stdlib.h>
#include "makespl.h"
#include "piv_ge_solver.h"

double H(int step , double x)
{
	if(step == 0)
		return 1;
	if(step == 1)
		return 2*x;
	return 2*x*H(step-1,x)-2*(step-1)*H(step-2,x);
}

double gdy_0(double a , int i , int k){}
double gdy_1(double a , int i , int k){}
double gdy_2(double a , int i , int k){}
double gdy_3(double a , int i , int k){}
double dla_y(double a){}

void make_spl(points_t *pts, spline_t *spl)
{
	

	int ILOSC = pts->n;
	matrix_t *eqs = NULL;
	double *x = pts->x;  //wartosci x-ow
	double *y = pts->y;  //wartosci y-ow
	int i , j , k;
	int m = 4;  //rozmier 
	double s = 0; //przechowywanie wzoru

	eqs = make_matrix(m,m+1);


	for(i = 0 ; i < m ; i++)
		{
				for(j = 0 ; j < m ; j++)
				{
					s = 0;
						for(k = 0 ; k < ILOSC ; k++)
						{
							if(j == 0)
								gdy_0(s , i , k);
							else if(j == 1)
								gdy_1(s , i , k);
							else if(j == 2)
								gdy_2(s , i , k);
							else
								gdy_3(s , i , k);
						}

					add_to_entry_matrix(eqs, i, j , s);
				}
			s=0;
			dla_y(s , ILOSC);

			add_to_entry_matrix(eqs , j , m , s);
		}

}				

