#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::generateBC1()
{
	PetscErrorCode ierr;
	PetscInt	mstart, nstart, m,n,i,j,M,N;
	PetscReal	**qx, **qy;
	PetscReal	**bc1x, **bc1y;
	Vec	bc1xGlobal, bc1yGlobal;
	PetscReal	nu=fluid.nu;
	PetscReal	alphaImplicit = 0.5;
	PetscReal	coeffMinus = 0.0, coeffPlus = 0.0;

	ierr = VecSet(bc1, 0.0); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(qPack, bc1, &bc1xGlobal, &bc1yGlobal); CHKERRQ(ierr);

	//U-FLUXES
	ierr = DMDAVecGetArray(uda, bc1xGlobal, &bc1x); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, qxLocal, &qx);	CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
	//X-faces
	coeffMinus = alphaImplicit*nu*2.0/dxU[0]/(dxU[0]+dxU[1]);
	coeffPlus = alphaImplicit*nu*2.0/dxU[M]/(dxU[M] + dxU[M-1]);

	for(j=nstart; j<nstart+n; j++)
	{
		//-X most nearby left edge, dirichlet b.c.
		if(mstart ==0) 	 bc1x[j][0] += coeffMinus*qx[j][-1]/fluid.dy; 
		
		//+X right edge, convective b.c. for ibm, dirichlet for cavity flow
		if(mstart+m-1 == M -1)
		{    
			bc1x[j][M-1] += coeffPlus*qx[j][M]/fluid.dy;
		}
	}

	//y-faces
	coeffMinus = alphaImplicit*nu*2.0/dyU[0]/(dyU[0]+dyU[1]);
	coeffPlus = alphaImplicit*nu*2.0/dyU[N]/(dyU[M] + dyU[M-1]);

	for(j=mstart; j<mstart+m; j++)
	{
		//-Y bottom edge, dirichlet b.c.
		if(nstart ==0)   bc1x[0][j] += coeffMinus*qx[-1][j];

		//+Y top edge, dirichlet b.c.
		if(nstart+n-1 == N -1)  bc1x[N-1][i] += coeffPlus * qx[M][i];
	}

	ierr = DMDAVecRestoreArray(uda, bc1xGlobal, &bc1x); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);

	//V-FLUXES
	ierr = DMDAVecGetArray(vda, bc1yGlobal, &bc1y); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);

	//x-faces
	coeffMinus = alphaImplicit*nu*2.0/dxV[0]/(dxV[0] + dxV[1]);
	coeffPlus = alphaImplicit*nu*2.0/dxV[M]/(dxV[M] + dxV[M-1]);

	for(j=nstart; j<nstart+n; j++)
	{
		//-X 
		if(mstart == 0)  bc1y[j][0] += coeffMinus*qy[j][-1];

		//+X
		if(mstart+m-1 == M-1) bc1y[j][M-1] += coeffPlus*qy[j][M]; 
	}
	//y-faces
	coeffMinus = alphaImplicit*nu*2.0/dyV[0]/(dyV[0] + dyV[1]);
	coeffPlus = alphaImplicit*nu*2.0/dyV[N]/(dyV[N] + dyV[N-1]);

	for(j=mstart; j<mstart+m; j++)
	{
		//-Y
		if(nstart == 0)	bc1y[0][j] += coeffMinus*qy[-1][j]/fluid.dx; 

		//+Y
		if(nstart+n-1 == N-1)  bc1y[N-1][j] += coeffPlus*qy[N][j]/fluid.dx; 
	}

	ierr = DMDAVecRestoreArray(vda, bc1yGlobal, &bc1y); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);
	ierr = DMCompositeRestoreAccess(qPack, bc1, &bc1xGlobal, &bc1yGlobal); CHKERRQ(ierr);

	return 0;
}
