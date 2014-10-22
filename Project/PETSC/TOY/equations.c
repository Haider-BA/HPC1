#include "headers.h"

#define  lap(m) ( (X[j][i+1].m - 2*X[j][i].m + X[j][i-1].m) / hx / hx + (X[j+1][i].m - 2*X[j][i].m + X[j-1][i].m) / hy / hy) 

PetscErrorCode equation(SNES snes, Vec soln_vec, Vec fun_vec, void* ptr)
{
	//data structure will be contained in ptr
	Data* data = (Data*) ptr;
	Params P = data->P;
	DA da = data->da;

	DALocalInfo info;
	PetscFunctionBegin;

	DAGetLocalInfo(da, &info);
	Vec sv_local; //need local vector containing ghost point information

	DACreateLocalVector(da, &sv_local);
	DAGlobalToLocalBegin(da, soln_vec, INSERT_VALUES, sv_local);
	DAGlobalToLocalEnd(da, soln_vec, INSERT_VALUES, sv_local);

	//to evaluate fuction, the sv_local soln_vec need to be converted into Field arrays
	Field **X;
	Field **F;

	DAVecGetArray(da, sv_local, &X);
	DAVecGetArray(da, fun_vec, &F);
	//since sv_local is local vector, X will have access to ghost points
	
	PetscInt i, j, Ny = P.Ny;
	PetscScalar hx = P.hx; hy = P.hy; x, y;

	for(i=info.xs; i<info.xs + info.xm; i++){
		for(j=info.ys; j<info.ys+info.ym; j++){
			x = i *hx;
			y = j *hy;

			if(j==0) {
				F[j][i].u = X[j][i].u - 0.5*(sin(2*PETSC_PI*x) + 1);
				F[j][i].v = X[j][i].v ;
			}
			else if (j== Ny-1){
				F[j][i].u = X[j][i].u;
				F[j][i].v = X[j][i].v - 0.5*(sin(2*PETSC_PI*x) + 1);
			}
			else {
				F[j][i].u = lap(u) + X[j][i].u * X[j][i].v ;
				F[j][i].v = lap(v) - X[j][i].u * X[j][i].v ;
			}
		}
	}

	DAVecRestoreArray(da, fun_vec, &F);
	DAVecRestoreArray(da, sv_local, &X);

	VecDestory(sv_local);

	PetscFuntionReturn(0);
}







