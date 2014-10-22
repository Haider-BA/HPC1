#include "petscsles.h"

int main(int argc, char** args)
{
	Vec x, b, u; //u exact solution
	Mat A;
	SLES sles;
	PetscRandom rctx; //random number generator context
	PetscReal norm;
	int i,j, I, J, Istart, Iend, ierr, m=4, n=4, its;
	PetscTruth flg;
	PetscScalar v,h, one = 1.0; neg_one = -1.0;
	KSP ksp;
	KSPType ksptype;
	PC pc;
	PCType pctype;
	PetscInitialize(&argc, &args, (char*)0, help);
	PetscOptionsGetInt(PETSC_NULL, "-m", &m, PETSC_NULL);
	PetscOptionsGetInt(PETSC_NULL, "-n", &n, PETSC_NULL);

	MatCreate(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, m*n, m*n, &A);
	MatSetFromOptions(A);
	MatGetOwnershipRange(A, &Istart, &Iend);
	for(I = Istart; I < Iend; I++){
		v = -1.0; i=I/n; j= I - i*n;
		if(i>0) { J = I-n; MatSetValues(A,1,&I, 1, &J, &v, INSERT_VALUES);}  //(I-1,J)
		if(i<m-1) {J=I+n; MatSetValues(A, 1, &I, 1, &J, &v, INSERT_VALUES);} //(I+1,J)
		if(j>0) { J = I-1; MatSetValues(A, 1, &I, 1, &J, &v, INSERT_VALUES);}//(I,J-1)		    
	       	if(j<n-1){J=I+1;  MatSetValues(A, 1, &I, 1, &J, &v, INSERT_VALUES);} //(I,J+1)
		v = 4.0; MatSetValues(A, 1, &I, 1, &I, &v, INSERT_VALUES); //(I,I)
	}
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

	VecCreate(PETSC_COMM_WORLD, &u);
	VecSetSizes(u, PETSC_DECIDE, m*n);
	VecSetFromOptions(u);
	VecDuplicate(u, &b);
	VecDuplicate(b, &x);

	//initializing values of b
	PetscOptionsHasName(PETSC_NULL, "-random_exact_sol", &flg);
	if(!flg){
		VecGetOwnershipRange(b, &Istart, &Iend);
		h = 1.0/(m+1);
		for(I=Istart; I<Iend; I++){
			v=0; i=I/n; j=I - i*n; h=1/(m+1);
			if(i==0) v = v + /*u(-1,j)*/ h*0;
			if(i==m-1) v = v + /* u(m,j) */ h *(m+1);
			if(j==0) v = v+ /*u(i,-1)*/ h*(i+1);
			if(j==n-1) v = v + /*u(i,n)*/ h*(i+1);

			if(v != 0) VecSetValues(b, 1, &I, &v, INSERT_VALUES);
			v = /* u(i,j) */ h*(i+1); VecSetValues(u, 1, &I, &v, INSERT_VALUES);
		}
		VecAssemblyBegin(b); VecAssemblyEnd(b);
		VecAssemblyBegin(u); VecAssemblyEnd(u);
		}else{
			PetscRandCreate(PETSC_COMM_WORLD, RANDOM_DEFAULT, &rctx);
			VecSetRandom(rctx, u); PetscRandomDestroy(rctx);
			MatMult(A,u,b);
		}
		PetscOptionsHasName(PETSC_NULL, "-view_exact_sol", &flg);
		if(flg) {VecView(u, PETSC_VIEWER_STDOUT_WORLD);}

	SLESCreate(PETSC_COMM_WORLD,&sles);
	SLESSetOperators(sles,A,A, DIFFERENT_NONZERO_PATTERN);
	SLESGetKSP(sles, &ksp);
	KSPSetTolerances(ksp, 1.e-2/((m+1)*(n+1)), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT);
	SLESSetFromOptions(sles);
	
	SLESSolve(sles, b, x, &its);

	PetscOptionsHasName(PETSC_NULL, "-view_sol_serial", &flg);
	if(flg) {VecView(x, PETSC_VIEWER_STDOUT_WORLD);}
	PetscOptionsHasName(PETSC_NULL, "-view_sol", &flg);
	if(flg){
		PetscScalar *xx;
		VectGetArray(x, &xx);
		VecGetOwnershipRange(x, &Istart, &Iend);
		PetscPrintf(PETSC_COMM_WORLD, "Solution Grid(without boundary conditions):\n");
		for(I=Istart; I<Iend; I++){
			i = I/n; j=I - i*n;
			PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%8.6f",xx[I-Istart]);
			if(j==(n-1)) PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n");
		}
		PetscSynchronizedFlush(PETSC_COMM_WORLD);
		VecRestoreArray(x, &xx);
	}

	#include "petscda.h"
	PetscOptionsHasName(PETSC_NULL, "-view_sol_x", &flg);
	if(flg) {
		PetscScalar *xx; DA da;
		AO ao; Vec x_da;
		DACreate2d(PETSC_COMM_WORLD, DA_NONPERIODIC, DA_STENCIL_STAR,n,m,PETSC_DECIDE, PETSC_DECIDE, 1,0, PETSC_NULL, PETSC_NULL, &da);
		DACreateGlobalVector(da, &x_da);
		DAGetAO(da, &ao);
		VecGetOwnershipRange(x, &Istart, &Iend);
		VecGetArray(x, &xx);
		for(I=Istart; I < Iend; I++){
			i=I; AOApplicationToPetsc(ao, 1, &i);
			VecSetValues(x_da, 1, &i, &&xx[I-Istart], INSERT_VALUES);
		}
		VecRestoreArray(x, &xx);
		VecAssemblyBegin(x_da); VecAssemblyEnd(x_da);

		PetscOptionsHasName(PETSC_NULL, "-view_sol_x_da", &flg);
		if(flg) VecView(x_da, PETSC_VIEWER_STDOUT_WORLD);
		VecView(x_da, PETSC_VIEWER_DRAW_(PETSC_COMM_WORLD));
		AODestroy(ao);
		DADestroy(da);
		VecDestroy(x_da);
	}

	//check solution
	VecAXPY(&neg_one, u, x);
	VecNorm(x, NORM_2, &norm);
	PetscPrintf(PETSC_COMM_WORLD, "Norm of error %A iterations %d\n", norm, its);
	SLESDestroy(sles);
	VecDestroy(u); VecDestroy(x);
	VecDestroy(b); MatDestroy(A);

	PetscFinalize();
}
