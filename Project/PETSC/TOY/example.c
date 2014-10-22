#include "headers.h"

PetscErrorCode guess(Vec guess, Data data);
PetscErrorCode equation(SNES snes, Vec SS, Vec FUN, void* ptr);
PetscErrorCode jacobian(SNES snes, Vec SS, Mat* jac, Mat* B, MatStructure* flag, void* ptr);

int main(int argc, char **argv)
{
	PetscInitalize(&argc, &argv, PETSC_NULL, PETSC_NULL);

	DA da;
	Params P;
	Data data;

	P = parameters(); //default construct 
	
	DACreate2d(PETSC_COMM_WORLD, DA_XPERIODIC, DA_STENCIL_STAR, P.Nx, P.Ny, PETSC_DECIDE, PETSC_DECIDE, 2,1,PETSC_NULL, PETSC_NULL, &da);
	data.P = P;
	data.da = da;

	Vec fun_vec, soln_vec; //preallocate space for function_vec, solution_vec
	Mat J;

	DACreateGlobalVector(da, &fun_vec);
	DACreateGlobalVector(da, &soln_vec);
	DAGetMatrix(da, MATMPIAIJ, &J);
// I need to say the details of these functions

	//MATMPIAIJ, the default way of storing sparse parallel matrices; and DAGetMatrix will determine location of nonzero elements in the Matrix
	
	SNES snes;
	SNESCreate(PETSC_COMM_WORLD, &snes);
	SNESSetFunction(snes, fun_vec, equation, &data);
	SNESSetJacobian(snes, J, J, jacobian, &data);


	KSP ksp;
	PC pc;
	SNESGetKSP(snes, &ksp);
	KSPSetType(ksp, KSPGMRES);
	KSPGMRESSetRestart(ksp, 100);
	KSPGetPC(KSP,&pc);
	PCSetType(pc, PCBJACOBI);
	PCSetFromOptions(pc);
	KSPSetFromOptions(ksp);
	SNESSetFromOptions(snes);

	//before solving, print information and seting guess initial
	guess(soln_vec, data);
	PetscPrintf(PETSC_COMM_WORLD, "grid size is %d x %d\n", P.Nx, P.Ny);
	//PetscPrintf will get root process to do print, if using printf, each process will do a print
	
	PetscLogDouble v1, v2, elapsed_time;
	PetscGetTime(&v1);
	SNESSolve(snes, PETSC_NULL, soln_vec);
	PetscGetTime(&v2);
	elapsed_time = v2 - v1;


	SNESConvergedReason snes_reas; //if SNES unable to find a solution, know what' wrong
	PetscInt lin_its, nonlin_its; // total number of GMRES iterations and number Newton iterations

	SNESGetConvergedReason(snes, &snes_reas);
	if(snes_res < 0){
		PetscPrintf(PETSC_COMM_WORLD, "SNES Failed! Reason %d\n", snes_reas);
	}

	SNESGetLinearSolveIterations(snes, &lin_its);
	SNESGetIterationNumber(snes, &nonlinear_its);
	PetscPrintf(pETSC_COMM_WORLD, "average number of GMRES iterations = (total GMRES its)/(Newton its) = %d %d = %.4e\n", lin_its, nonlin_its, ((double) lin_its)/nonlin_its);

	PetscPrintf(PETSC_COMM_WORLD, "nonlinear eqns sovled in %.2f s\n", elapsed_time);

	FILE *fid;
	printvec(soln_vec, "soln.dat");

	PetscFOpen(PETSC_COMM_WORLD, "soln_info.dat", "w", &fid);
	PetscFPrintf(PETSC_COMM_WORLD, fid, "%d\n%d", P.Nx, p.Ny);
	PetscFClose(PETSC_COMM_WORLD, fid);

	VecDestroy(soln_vec);
	VecDestroy(fun_vec);
	MatDestroy(J);

	PetscFinalize();
	return 0;
}


