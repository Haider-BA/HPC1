#include "../NavierStokesSolver.h"

PetscErrorCode	NavierStokesSolver::generateBC2()
{
	PetscErrorCode ierr;
	PetscInt	mstart, nstart, m, n, i, j;
	PetscReal	**qx, **qy;
	PetscReal	**bc22;
	Vec	bc22Global;

	ierr = VecSet(bc2, 0.0); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(lambdaPack, bc2, &bc22Global); CHKERRQ(ierr);
	
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(pda, bc22Global, &bc22); CHKERRQ(ierr);
	ierr = DMDAGetCorners(pda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);

	//x-faces
	for(j=nstart; j<nstart+n; j++)
	{
		//-X
		if(mstart ==0)	bc22[j][0] -= qx[j][-1];
		//+X
		if(mstart+m-1 == fluid.nx-1) bc22[j][fluid.nx-1] += qx[j][fluid.nx-1];
	}

	//y-faces
	for(j=mstart; j<mstart+m; j++)
	{
		//-Y
		if(nstart ==0) bc22[0][j] -= qy[-1][j];
		//+Y
		if(nstart+n-1 == fluid.ny-1)	bc22[fluid.nx-1][j] += qy[fluid.nx-1][j];
	}
	
	ierr = DMDAVecRestoreArray(pda, bc22Global, &bc22); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(lambdaPack, bc2, &bc22Global); CHKERRQ(ierr);

	return 0;
}
