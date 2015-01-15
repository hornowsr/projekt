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
 
void sprawdz(spline_t *spl, matrix_t *eqs)
{	
	int m=4;
	if(alloc_spl(spl,1)==0)
	{
		double a0=get_entry_matrix(eqs, 0 , m);
		double a1=get_entry_matrix(eqs, 1 , m);
		double a2=get_entry_matrix(eqs, 2 , m);
		double a3=get_entry_matrix(eqs, 3 , m);
		printf("%lf, %lf, %lf, %lf\n",a0,a1,a2,a3);
		double x;

		spl->x[0]=x;
		spl->f[0]=a0+a1*H(1,x)+a2*H(2,x)+a3*H(3,x);
		spl->f1[0]=2*a1+a2*8*x*a3*(24*x*x-12);
		spl->f2[0]=a2*8+a3*48*x;
		spl->f3[0]=a3*48;
	}
}
double dla_y( int j , int max, double *z, double *c)
{
	int k;
	double s=0;
	for(k=0 ; k<max ; k++)
	{
		if(j == 0)
			s = s+c[k];
		if(j == 1)
			s = s+c[k]*H(1,z[k]);
		if(j == 2)
			s = s+c[k]*H(2,z[k]);
		if(j == 3)
			s = s+c[k]*H(3,z[k]);
	}
	return s;
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
	if(j == 0)
		s = s+H(1,x[k]);
	else if(j == 1)
		s = s+H(1,x[k])*H(1,x[k]);
	else if(j == 2)
		s = s+H(2,x[k])*H(1,x[k]);
	else if(j == 3)
		s = s+H(3,x[k])*H(1,x[k]);

	return s;
}

double gdy_2(double s, int j, int k, double *x)
{
	if(j == 0)
		s = s+H(2,x[k]);
	else if(j == 1)
		s = s+H(1,x[k])*H(2,x[k]);
	else if(j == 2)
		s = s+H(2,x[k])*H(2,x[k]);
	else if(j == 3)
		s = s+H(3,x[k])*H(2,x[k]);

	return s;
}

double gdy_3(double s, int j , int k , double *x)
{
	if(j == 0)
		s = s+H(3,x[k]);
	else if(j == 1)
		s = s+H(1,x[k])*H(3,x[k]);
	else if(j == 2)
		s = s+H(2,x[k])*H(3,x[k]);
	else if (j == 3)
		s = s+H(3,x[k])*H(3,x[k]);
	
	return s;
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

	sprawdz(spl,eqs);
}

