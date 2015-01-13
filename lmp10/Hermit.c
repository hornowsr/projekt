#include <stdio.h>
#include <stdlib.h>


double H(int step , double x)
{
	if(step == 0)
		return 1;
	if(step == 1)
		return 2*x;
	return 2*x*H(step-1,x)-2*(step-1)*H(step-2,x);
}

