#include "../TcSolver.h"

PetscErrorCode	TcSolver::calculateForce()
{
	PetscErrorCode ierr;
	PetscInt	mstart, nstart, m, n;
	Vec	fGlobal, fxGlobal, fyGlobal;
	PetscReal	**fx, **fy, forceOnProcess[2];

	ierr = DMCompositeGetAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);
	ierr = MatMult(ET, fGlobal, regularizedForce); CHKERRQ(ierr);

	ierr = DMCompositeGetAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal); CHKERRQ(ierr);

	//x-direction
	ierr = DMDAVecGetArray(uda, fxGlobal, &fx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	forceOnProcess[0] = 0;
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			forceOnProcess[0] += fluid.dy * fx[j][i];
		}
	}
	ierr = DMDAVecRestoreArray(uda,fxGlobal,&fx); CHKERRQ(ierr);

	//y-direction
	ierr = DMDAVecGetArray(vda, fyGlobal, &fy); CHKERRQ(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n,NULL); CHKERRQ(ierr);
	forceOnProcess[1] = 0;
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			forceOnProcess[1] += fluid.dy*fy[j][i];
		}
	}
	ierr = DMDAVecRestoreArray(vda, fyGlobal, &fy); CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(qPack, regularizedForce, &fxGlobal, &fyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);

	ierr = MPI_Reduce(forceOnProcess, force, 2, MPIU_REAL, MPI_SUM, 0, MPI_COMM_WORLD); CHKERRQ(ierr);

	return 0;
}
