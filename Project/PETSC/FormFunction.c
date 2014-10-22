#include"petscda.h"
#include"petscsnes.h"
#include<math.h>

int FormBratuFunction(SNES snes, Vec v, Vec f, void* ctx)
{
	UserBratuCtx* bratu = (UserBratuCtx *) ctx;
	DA da = bratu->da;
	double lambda = bratu->lambda;
	double h = bratu->h;
	Vec lv;
	int i,j;
	int lli, llj, ni, nj;
	const double **varr;
	double **fvarr;

	DAGetCorners(da, &lli, &llj, 0, &ni, &nj, 0);
	DAGetLocalVector(da, &lv);
	DAGlobalToLocalBegin(da, v, INSERT_VALUES, lv);
	DAGlobalToLocalEnd(da, v, INSERT_VALUES,lv);

	DAVecGetArray(da, lv, (void**) &varr);
	DAVecGetArray(da, f, (void**) &fvarr);

	for(j=llj; j<llj+nj; j++){
		for(i=lli; i<lli+ni; i++){
			if(i==0 || j==0 ||i=bratu->n+1 || j==bratu->n+1){
				fvarr[j][i] = 0.0;
			}
			else{
				fvarr[j][i] = -(varr[j-1][i] + varr[j][i-1] + varr[j+1][i] + varr[j][i+1] - 4*varr[j][i])/(h*h) - lambda*exp(varr[j][i]);
			}
		}
	}
	DAVecRestoreArray(da, f, (void**) &fvarr);
	DAVecRestoreArray(da, lv, (void**) &varr);
	DARestoreLocalVector(da, &lv);

	return 0;
}
