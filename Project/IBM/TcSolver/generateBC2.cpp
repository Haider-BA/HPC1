#include "../TcSolver.h"

	PetscErrorCode TcSolver::generateBC2()
{
	PetscErrorCode	ierr;
	PetscInt	mstart, nstart, m, n, i, j;
	PetscReal	**qx, **qy;
	PetscReal	**bc22;
	Vec		bc2Global;

	ierr = VecSet(bc2, 0.0); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(lambdaPack, bc2, &bc2Global, NULL); CHKERRQ(ierr);
	
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);

	ierr = DMDAVecGetArray(pda, bc2Global, &bc22); CHKERRQ(ierr);
	ierr = DMDAGetCorners(pda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);

	//x-faces
	
	for(j=nstart; j<nstart+n; j++)
	{
			//-X
		if(mstart == 0) //if -X is current process
		{
			bc22[j][0] -= qx[j][-1];
		}
		//+X
		if(mstart+m-1 == fluid.nx-1)
		{
			bc22[j][fluid.nx-1] += qx[j][fluid.nx-1];
		}
	}
	

	//y-faces
	for(i=mstart; i<mstart+m; i++)
	{
		//-Y
		if(nstart ==0)
		{
			bc22[0][i] -= qy[-1][i];
		}
		//+Y
		if(nstart+n-1 == fluid.ny-1)
		{
			bc22[fluid.ny-1][i] += qy[fluid.ny-1][i];
		}
	}
	
	ierr = DMDAVecRestoreArray(pda, bc2Global, &bc22); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(lambdaPack, bc2, &bc2Global, NULL); CHKERRQ(ierr);

	return 0;
}





