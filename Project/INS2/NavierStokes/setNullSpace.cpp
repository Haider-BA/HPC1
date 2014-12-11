/*
 * the null space for Poisson system in a flow with no immersed boundary is the constant vector of size equal to the number of pressure variables
 *
 * such a null space can be specified and automaticaly hanled by PTESc by passing 
 * PETSC_TRUE as the second parameter in MatNullSpaceCreate
 */
#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::setNullSpace()
{
	PetscErrorCode ierr;
	MatNullSpace nsp;

	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, NULL, &nsp); CHKERRQ(ierr);
	ierr = KSPSetNullSpace(ksp2, nsp);
	ierr = MatNullSpaceDestroy(&nsp);

	return 0;
}
