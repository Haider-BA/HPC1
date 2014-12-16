#include "../TcSolver.h"

PetscErrorCode TcSolver::initializeLambda()
{
	PetscErrorCode ierr;
	Vec phi, fTilde;

	ierr = DMCompositeGetAccess(NavierStokesSolver::lambdaPack, NavierStokesSolver::lambda, &phi, &fTilde);	CHKERRQ(ierr);


	ierr = DMCompositeRestoreAccess(NavierStokesSolver::lambdaPack, NavierStokesSolver::lambda, &phi, &fTilde);

	return 0;
}
