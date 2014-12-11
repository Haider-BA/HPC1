
#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::initializeLambda()
{
	PetscErrorCode ierr;
	Vec	phi;

	ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi);

	VecZeroEntries(phi);

	return 0;
}
