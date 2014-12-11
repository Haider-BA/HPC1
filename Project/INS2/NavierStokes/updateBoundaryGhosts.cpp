#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::updateBoundaryGhosts()
{
	PetscErrorCode ierr;
	PetscInt	mstart, nstart, m, n, i, j, M, N;
	PetscReal	**qx, **qy;
	PetscReal	dt = fluid.dt;

	//u-fluxes
	ierr = DMDAVecGetArray(uda, qxLocal, &qx);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL);
	ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
		
	for(j=nstart; j<nstart+n; j++)
	{
		//left-edge(inlet)
		if(mstart==0) //global indices
		{		//dirichlet b.c.
			qx[j][-1] = fluid.bc[0][0] * fluid.dy;
		}

		//right-edge(outlet)
		if(mstart+m == M) // M global indices range in x direction
		{ //convective approximation
			PetscReal beta = fluid.bc[0][3] * dt / dxU[M];
			qx[j][M] = (1-beta)*qx[j][M] + beta*qx[j][M-1];
		}	
	}
	
	for(i=mstart; i<mstart+m; i++)
	{
		//bottom wall
		if(nstart==0)
		{ //dirichlet b.c.
			qx[-1][i] = fluid.bc[0][1] * fluid.dx; //take care 
		}

		//top wall
		if(nstart+n ==N)
		{ //dirichlet b.c.
			qx[N][i] = fluid.bc[0][2] * fluid.dx;
		}
	}
	
	ierr = DMDAVecRestoreArray(uda, qxLocal, &qx);

	//v-fluxes
	ierr = DMDAVecGetArray(vda, qyLocal, &qy);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL);
	ierr = DMDAGetInfo(vda, NULL, &M, &N, NULL, NULL, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
		
	for(j=nstart; j<nstart+n; j++)
	{
		//left-edge(inlet)
		if(mstart==0) //global indices
		{		//dirichlet b.c.
			qy[j][-1] = fluid.bc[1][0] * fluid.dy;
		}

		//right-edge(outlet)
		if(mstart+m == M) // M global indices range in x direction
		{ //convective approximation
			PetscReal beta = fluid.bc[1][3] * dt / dxV[M];
			qy[j][M] = (1-beta)*qy[j][M] + beta*qy[j][M-1] ;
		}	
	}
	
	for(i=mstart; i<mstart+m; i++)
	{
		//bottom wall
		if(nstart==0)
		{ //dirichlet b.c.
			qy[-1][i] = fluid.bc[1][1] * fluid.dx; //take care, why GWU didn't multiply dx here 
		}

		//top wall
		if(nstart+n ==N)
		{ //dirichlet b.c.
			qy[N][i] = fluid.bc[1][2] * fluid.dx;
		}
	}
	ierr = DMDAVecRestoreArray(vda, qyLocal, &qy);
	
	return 0;
}


