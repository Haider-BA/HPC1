/*
 * set default options for Kyrlov solvers
 *
 * ksp1, used solving intermediate velocity, default type CG
 * relative tolerance for convergence is 10^-5
 * initial guess for the solution is obtained from the output vector supplied
 * command line arguments to set options for this solver must have the prefix 'sys1_'
 *
 * ksp2, used solving pressure, same set as ksp1
 */

#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::createKSPs()
{
	PetscErrorCode ierr;

	//linear system for intermmediate velocity
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp1); CHKERRQ(ierr);
//	ierr = KSPSetOptionPrefix(ksp1, "sys1_"); CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp1, 1e-5, 0.0, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp1, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp1, PETSC_TRUE); CHKERRQ(ierr);
	ierr = KSPSetType(ksp1, KSPCG); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp1); CHKERRQ(ierr);

	//linear system for Poisson solver
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp2); CHKERRQ(ierr);
//	ierr = KSPSetOptionPrefix(ksp2, "sys2_"); CHKERRQ(ierr);
	ierr = KSPSetTolerances(ksp1, 1e-5, 0.0, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp2, QTBNQ, QTBNQ, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
	ierr = KSPSetInitialGuessNonzero(ksp2, PETSC_TRUE); CHKERRQ(ierr);
	ierr = KSPSetType(ksp2, KSPCG); CHKERRQ(ierr);
	ierr = KSPSetFromOptions(ksp2); CHKERRQ(ierr);

	return 0;
}
	
