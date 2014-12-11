#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver :: createLocalToGlobalMappingsLambda()
{
	PetscErrorCode ierr;
	PetscInt	m,n,i,j,mstart,nstart;
	PetscReal	**lp;
	PetscInt	localIdx;

	ierr = VecGetOwnershipRange(lambda, &localIdx, NULL); CHKERRQ(ierr);
//populate local vector with global indices
//values outside the domain are never accessed and not set
//P
	ierr = DMCreateLocalVector(pda, &pMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(pda, pMapping, &lp); CHKERRQ(ierr);
	ierr = DMDAGetCorners(pda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i < mstart+m; i++)
		{
			lp[j][i] = localIdx;
			localIdx++;
		}
	}

	ierr = DMDAVecRestoreArray(pda, pMapping, &lp); CHKERRQ(ierr);

	ierr = DMDALocalToLocalBegin(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);
	ierr = DMDALocalToLocalEnd(pda, pMapping, INSERT_VALUES, pMapping); CHKERRQ(ierr);
	return 0;
}

