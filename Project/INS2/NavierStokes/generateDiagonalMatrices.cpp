#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver :: generateDiagonalMatrices()
{
	PetscErrorCode ierr;
	PetscInt	mstart, nstart, m,n,i,j;
	Vec	MHatxGlobal, MHatyGlobal;
	Vec	RInvxGlobal, RInvyGlobal;
	Vec	BNxGlobal, BNyGlobal;
	PetscReal	**MHatx, **MHaty;
	PetscReal	**RInvx, **RInvy;
	PetscReal	**BNx, **BNy;

	ierr = DMCompositeGetAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(qPack, BN, &BNxGlobal, &BNyGlobal); CHKERRQ(ierr);

	//x-direction
	ierr = DMDAVecGetArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, BNxGlobal, &BNx);  CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i <mstart+m; i++)
		{
//			MHatx[j][i] = (i<fluid.nx-1)?0.5*(mesh->dx[i] + mesh->dx[i+1]):0.5*(mesh->dx[i] + mesh->dx[0]);
			MHatx[j][i] = fluid.dx; 
			RInvx[j][i] = 1.0/fluid.dy;
			BNx[j][i] = fluid.dt/(MHatx[j][i]*RInvx[j][i]);
		}
	}

	ierr = DMDAVecRestoreArray(uda, MHatxGlobal, &MHatx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uda, RInvxGlobal, &RInvx); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uda, BNxGlobal, &BNx);   CHKERRQ(ierr);

	//y-direction
	ierr = DMDAVecGetArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, BNyGlobal, &BNy); CHKERRQ(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
//			MHaty[j][i] = (j<mesh->ny-1) ? 0.5*(mesh->dy[j] + mesh->dy[j+1]):0.5*(mesh->dy[j]+mesh->dy[0]);
			MHaty[j][i] = fluid.dy;
			RInvy[j][i] = 1.0/fluid.dx;
			BNy[j][i] = fluid.dt/(MHaty[j][i]*RInvy[j][i]);
		}
	}

	ierr = DMDAVecRestoreArray(vda, MHatyGlobal, &MHaty); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(vda, RInvyGlobal, &RInvy); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(vda, BNyGlobal, &BNy); CHKERRQ(ierr);

	ierr = DMCompositeRestoreAccess(qPack, MHat, &MHatxGlobal, &MHatyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(qPack, RInv, &RInvxGlobal, &RInvyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(qPack, BN, &BNxGlobal, &BNyGlobal); CHKERRQ(ierr);

	return 0;
}

