/*
 * vectors stored as distributed arrays can be accessed using multi-dimensional arrays on every process
 * with each index referring to the numbering along each Cartesian direction
 * the elements of the vector also have a global ordering
 *
 * this function generate the map from multi-dimensional indexing of each of the local flux vectors "qx, qy" to the global indices of the composite flux vector "q"
 *
 */
#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::createLocalToGlobalMappingsFluxes()
{
	PetscErrorCode ierr;
	PetscInt m, n, i, j, mstart, nstart;
	PetscReal	**lx, **ly;
	PetscInt	localIdx;

	ierr = VecGetOwnershipRange(q, &localIdx, NULL); CHKERRQ(ierr);
//populate local vectors with global indices
//set values to -1 if the cell is outside the domain
//U
	ierr = DMCreateLocalVector(uda, &uMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(uda, uMapping, &lx); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m;i++)
		{
			lx[j][i] = -1;
			if(i>mstart && i<mstart+m-1 && j>nstart && j<nstart+n-1)
			{
				lx[j][i] = localIdx;
				localIdx++;
			}
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &lx); CHKERRQ(ierr);

//V
	ierr = DMCreateLocalVector(vda, &vMapping); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(vda, vMapping, &ly); CHKERRQ(ierr);
	ierr = DMDAGetGhostCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m;i++)
		{
			ly[j][i] = -1;
			if(i>mstart && i<mstart+m-1 && j>nstart && j<nstart+n-1)
			{
				ly[j][i] = localIdx;
				localIdx++;
			}
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &ly); CHKERRQ(ierr);

//scatter from local to local to obtain correct values in ghost cells

ierr = DMDALocalToLocalBegin(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);
ierr = DMDALocalToLocalEnd(uda, uMapping, INSERT_VALUES, uMapping); CHKERRQ(ierr);

ierr = DMDALocalToLocalBegin(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);
ierr = DMDALocalToLocalEnd(vda, vMapping, INSERT_VALUES, vMapping); CHKERRQ(ierr);

return 0;
}



