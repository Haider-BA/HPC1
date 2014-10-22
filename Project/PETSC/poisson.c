#include<math.h>
#include "petscksp.h"

extern Mat FormLaplacian2d(int);
extern Vec FormVecFromFunction2d(int, double (*)(double, double));

double func(double x, double y){
	return sin(x*M_PI)*sin(y*M_PI);
}

int main(int argc, char* argv[])
{
	KSP sles;
	Mat A;
	Vec b, x;
	int its, n;

	PetscInitialize(&argc, &argv, 0, 0);
	n = 10;
	PetscOptionsGetInt(PETSC_NULL, "-n", &n, 0);

	A = FormLaplacian2d(n);
	b = FormVecFromFunction2d(n, func);
	VecDuplicat(b, &x);
	KSPCreate(PETSC_COMM_WORLD, &sles);
	KSPSetOperators(sles,A,A,DIFFERENT_NONZERO_PATTERN);
	KSPSetFromOptions(sles);
	KSPSolve(sles, b, x);
	KSPGetInterationNumber(sles, &its);
	PetscPrintf(PETSC_COMM_WORLD, "solution in %d iterations is: \n");
	VecView(x, PETSC_VIEWER_STDOUT_WORLD);

	MatDestroy(A); VecDestory(b); VecDestroy(x);
	KSPDestroy(sles);
	PetscFinalize();
	return 0;
}

