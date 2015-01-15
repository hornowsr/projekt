#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

double H( int stopien, double x )
{
        if( stopien == 0 ) return 1;
        if( stopien == 1 ) return 2*x;
        return 2*x * H( stopien-1, x) - 2*(stopien-1)* H(stopien-2,x);
}
 
double dla_y( int j , int max, double *z, double *c)
{
}

double gdy_0(double s, int j, int k , double *x)
{
	if(j == 0)
		s = s+H(0,1);
	else if(j == 1)
		s = s+H(1,x[k]);
	else if(j == 2)
		s = s+H(2,x[k]);
	else if(j == 3)
		s = s+H(3,x[k]);

	return s;
}

double gdy_1(double s, int j , int k , double *x)
{
}

double gdy_2(double s, int j, int k, double *x)
{
}

double gdy_3(double s, int j , int k , double *x)
{
}
void
make_spl(points_t * pts, spline_t *spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	int		i, j, k, max=pts->n;
	int 		m = 4;
	double 		s=0;
	eqs = make_matrix(m, m+1);
	
	for (j = 0; j < m; j++) 
	{
		
		for (i = 0; i < m; i++)
		{
			s=0;
			for (k = 0; k < max; k++)
			{	
				if(i==0)
				s=gdy_0(s,j,k,x);
				else if (i==1)
				s=gdy_1(s,j,k,x);
				else if (i==2)
				s=gdy_2(s,j,k,x);
				else if (i==3)
				s=gdy_3(s,j,k,x);
			}
			add_to_entry_matrix(eqs, j ,i,s);
		}
		s=0;
		s=dla_y(j,max,x,y);
			add_to_entry_matrix(eqs, j, m,s );
	}
	if(piv_ge_solver(eqs))
	{
		spl->n = 0;
		return;
	}

}

