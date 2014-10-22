#include"petsc.h"
#include"petscvec.h"
#include"petscda.h"

vec FormVecFromFunctionDA2d(DA grid, int n, double(*f) (double,double) )
{
	Vec V;
	int is, ie, js, je, in, jn, i, j;
	double h;
	double **vval;

	h = 1.0/(n+1);
	DACreateGlobalVector(grid, &V);
	DAVecGetArray(grid, V, (void**) &vval);

	DAGetCorners(grid, &is, &js, 0, &in, &jn, 0);
	ie = is + in - 1;
	je = js + jn - 1;
	for(i=is; i<=ie; i++){
		for(j=js; j<=je; j++)
		{
			vval[j][i] = (*f) ( (i+1)*h, (j+1)*h);
		}
	}
	DAVecRestoreArray(grid, V, (void**) &vval);

	return V;
}

