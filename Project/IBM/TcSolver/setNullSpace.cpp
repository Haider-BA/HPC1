#include "../TcSolver.h"

PetscErrorCode	TcSolver::setNullSpace()
{
	PetscErrorCode ierr;
	MatNullSpace	nsp;
	Vec	phiPortion;
	PetscInt	numPhi;

	ierr = VecSet(nullSpaceVec, 0.0);
	ierr = DMCompositeGetAccess(NavierStokesSolver::lambdaPack, nullSpaceVec, &phiPortion, NULL);
	ierr = VecGetSize(phiPortion, &numPhi);
	ierr = VecSet(phiPortion, 1.0/sqrt(numPhi));
	ierr = DMCompositeRestoreAccess(NavierStokesSolver::lambdaPack,nullSpaceVec,&phiPortion,NULL);
	
	ierr = MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,1,&nullSpaceVec,&nsp);
	ierr = KSPSetNullSpace(NavierStokesSolver::ksp2,nsp);
	ierr = MatNullSpaceDestroy(&nsp);

	return 0;
}
