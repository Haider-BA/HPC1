inline PetscReal du2dx2(PetscReal uMinus, PetscReal uCenter, PetscReal uPlus, PetscReal dxMinus, PetscReal dxPlus)
{
	return (dxPlus*uMinus + dxMinus*uPlus - (dxPlus + dxMinus)*uCenter)*2.0/dxMinus/dxPlus/(dxMinus+dxPlus);
}

 #include "../NavierStokesSolver.h"

/*
 * calculate explicit terms in the discretized Navier-Stokes equations
 *
 * include convection term and explicit portion of the diffusion term, the velocity value at the previous time step that appears in the time discretization is also added to the explicit terms
 *
 */

PetscErrorCode	NavierStokesSolver::calculateExplicitTerms()
{
	PetscErrorCode ierr;
	PetscInt	mstart, nstart, m,n,i,j,M,N;
	Vec		HxGlobal, HyGlobal;
	Vec		rxGlobal, ryGlobal;
	PetscReal	**qx, **qy;
	PetscReal	**Hx, **Hy, **rx, **ry;
	PetscReal	HnMinus1, u, v;
	PetscReal	uNorth, uEast, uWest, uSouth;
	PetscReal	vNorth, vEast, vWest, vSouth;
	PetscReal	convectionTerm, diffusionTerm;
	PetscReal	nu = fluid.nu;
	PetscReal	alphaExplicit = 0.5,
			gamma = 1.5,
			zeta = -0.5;
	PetscReal	dt= fluid.dt;

	//copy fluxes to local vectors
	ierr = DMCompositeScatter(qPack, q, qxLocal, qyLocal); CHKERRQ(ierr);
	
	ierr = DMCompositeGetAccess(qPack, H, &HxGlobal, &HyGlobal); CHKERRQ(ierr);
	ierr = DMCompositeGetAccess(qPack, rn, &rxGlobal, &ryGlobal); CHKERRQ(ierr);

	//access local vectors through multi-dimensional pointers
	ierr = DMDAVecGetArray(uda, qxLocal, &qx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, qyLocal, &qy); CHKERRQ(ierr);

	//x-component
	ierr = DMDAVecGetArray(uda, HxGlobal, &Hx); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, rxGlobal, &rx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	ierr = DMDAGetInfo(uda, NULL, &M, &N, NULL, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m;i++)
		{
			//convection term
			u = qx[j][i]/fluid.dy;
			uWest = 0.5*(u+qx[j][i-1]/fluid.dy);
			uEast = 0.5*(u+qx[j][i+1]/fluid.dy);

			//first check if the node is adjacent to -Y / +Y
			uSouth = (j>0) ? 0.5*(u+qx[j-1][i]/fluid.dy):qx[j-1][i]/fluid.dx;
			uNorth = (j<N-1)?0.5*(u+qx[j+1][i]/fluid.dy):qx[j+1][i]/fluid.dx;

			vSouth = 0.5*(qy[j-1][i]/dxU[i] + qy[j-1][i+1]/dxU[i+1]);
			vNorth = 0.5*(qy[j][i]/dxU[i] + qy[j][i+1]/dxU[i+1]);

			//Hx = d(u^2)/dx + d(uv)/dy
			HnMinus1 = Hx[j][i];
			Hx[j][i] = (uEast*uEast - uWest*uWest)/(0.5*(dxU[i]+ dxU[i+1])) + (uNorth*vNorth - uSouth*vSouth)/fluid.dy;
			convectionTerm = gamma*Hx[j][i] + zeta*HnMinus1;

			//diffusion term
			uWest = qx[j][i-1]/fluid.dy;
			uEast = qx[j][i+1]/fluid.dy;
			uSouth = (j>0) ? qx[j-1][i]/fluid.dy: qx[j-1][i];
			uNorth = (j<N-1)?qx[j+1][i]/fluid.dy: qx[j+1][i];

			//Dx = d^2(u)/dx^2 + d^2(u)/dy^2
			diffusionTerm = alphaExplicit*nu*(du2dx2(uWest,u,uEast, dxU[i], dxU[i+1]) + du2dx2(uSouth, u, uNorth, dyU[j], dyU[j+1]));

			rx[j][i] = (u/dt - convectionTerm + diffusionTerm);
		}
	}		
		ierr = DMDAVecRestoreArray(uda, HxGlobal, &Hx); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(uda, rxGlobal, &rx); CHKERRQ(ierr);

// y-component

		ierr = DMDAVecGetArray(vda, HyGlobal, &Hy); CHKERRQ(ierr);
		ierr = DMDAVecGetArray(vda, ryGlobal, &ry); CHKERRQ(ierr);
		ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
		for(j=nstart; j<nstart+n; j++)
		{
			for(i=mstart; i<mstart+m; i++)
			{
				//convection term
				v = qy[j][i]/fluid.dx;
				vSouth = 0.5*(v + qy[j-1][i]/fluid.dx);
				vNorth = 0.5*(v + qy[j+1][i]/fluid.dx);

				uWest = 0.5*(qx[j][i-1]/dyV[j] + qx[j+1][i-1]/dyV[j+1]);
				uEast = 0.5*(qx[j][i]/dyV[j] + qx[j+1][i]/dyV[j+1]);
				vWest = (i>0) ? 0.5*(v + qy[j][i-1]/fluid.dx):qy[j][i-1];
				vEast = (i<M-1)?0.5*(v+qy[j][i+1]/fluid.dx):qy[j][i+1];

				//Hx = d(uv)/dx + d(v^2)/dy
				HnMinus1 = Hy[j][i];
				Hy[j][i] = (uEast*vEast - uWest*vWest)/fluid.dx
					+ (vNorth*vNorth - vSouth*vSouth)/(0.5*(dyV[j] + dyV[j+1]));
				convectionTerm = gamma * Hy[j][i] + zeta * HnMinus1;
				//diffusion term
				vSouth = qy[j-1][i]/fluid.dx;
				vNorth = qy[j+1][i]/fluid.dx;
				vWest = (i>0)? qy[j][i-1]/fluid.dx: qy[j][i-1];
				vEast = (i<M-1) ? qy[j][i+1]/fluid.dx:qy[j][i+1];

				//Dy = d^2(v)/dx^2 + d^2(v)/dy^2
				diffusionTerm = alphaExplicit*nu*(du2dx2(vWest, v, vEast, dxV[i], dxV[i+1])+ du2dx2(vSouth, v, vNorth, dyV[j], dyV[j+1]));

				ry[j][i] = (v/dt - convectionTerm + diffusionTerm);
			}
		}

		ierr = DMDAVecRestoreArray(vda, HyGlobal, &Hy); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(vda, ryGlobal, &ry); CHKERRQ(ierr);
		
		ierr = DMDAVecRestoreArray(uda, qxLocal, &qx); CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(vda, qyLocal, &qy); CHKERRQ(ierr);

		ierr = DMCompositeRestoreAccess(qPack, H, &HxGlobal, &HyGlobal); CHKERRQ(ierr);
		ierr = DMCompositeRestoreAccess(qPack, rn, &rxGlobal, &ryGlobal); CHKERRQ(ierr);

		return 0;
	}


