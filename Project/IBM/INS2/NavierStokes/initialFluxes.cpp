#include "../NavierStokesSolver.h"

PetscErrorCode	NavierStokesSolver::initializeFluxes()
{
	PetscErrorCode ierr;
	PetscInt mstart, nstart, m, n;
	PetscReal	**qx, **qy;
	Vec	qxGlobal, qyGlobal;

	//read in initVel
	PetscReal	initVel[2] = {fluid.initialVelocity[0], fluid.initialVelocity[1]};

	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal);
	//u-fluxes
	ierr = DMDAVecGetArray(uda, qxGlobal, &qx);
	ierr = DMDAGetCorners(uda, &mstart, &nstart,NULL,&m,&n,NULL);
	for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
				{
					qx[j][i] = initVel[0]*fluid.dy;
				}
		}
	ierr = DMDAVecGetArray(uda,qxGlobal,&qx);

	//v-fluxes
	ierr = DMDAVecGetArray(vda, qyGlobal, &qy);
	ierr = DMDAGetCorners(vda, &mstart, &nstart,NULL,&m,&n,NULL);
	for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
				{
					qy[j][i] = initVel[1]*fluid.dx;
				}
		}
	ierr = DMDAVecGetArray(vda,qyGlobal,&qy);

	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal);

	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal);

	return 0;
	}


