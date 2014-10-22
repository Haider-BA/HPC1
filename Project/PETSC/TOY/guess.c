#include "headers.h"

PetscErrorCode guess(Vec guess, Data data)
{
	PetscFunctionBegin;
//for auxiliary functions that do not call PetscInitialize, should do this first
	DA da = data.da;
	DALocalInfo info; //to obtain local grid information
	Field **g;

	DAGetLocalInfo(da, &info);
	DAVecGetArray(da, guess, &g); //guess is a global DA vector

	PetscInt i,j;
	PetscScalar x, y;

	for(i=info.xs; i<info.xs + info.xm; i++){
		for(j=infor.ys; j<info.ys+info.ym; j++){

			x = i*hx;
			y = j*hy;

			g[j][i].u = 0.5*(sin(2*PETSC_PI*x) + 1) * pow(1-y,2);
			g[j][i].v = 0.5*(sin(2*PETSC_PI*x) + 1) * pow(y,2);
		}
	}

	DAVecRestoreArray(da, guess, &g);

	PetscFunctionReturn(0);  //for error handling, if else using C "return"

}
