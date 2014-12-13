 /* the size of array of each of the flow variables are:
 *
 * sizeof(u)= (nx-1) * ny
 * sizeof(v)= nx * (ny-1)
 * sizeof(p) = nx * ny
 *
*/

#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::createDMs()
{
	PetscErrorCode ierr;
	PetscInt 	m,n;
	const PetscInt *lxp, *lyp;
	PetscInt	*lxu, *lyu, *lxv, *lyv;
	PetscInt	numX, numY;

	//Create distributed array data structures
	// pressure
	numX = fluid.nx;
	numY = fluid.ny;
	ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_GHOSTED, DMDA_STENCIL_STAR, numX, numY, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &pda); CHKERRQ(ierr);
	ierr = DMDAGetOwnershipRanges(pda, &lxp, &lyp, NULL); CHKERRQ(ierr);
//obtain local lenght of p in x, y direction
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
//obtain num of procs in x, y direction and later mapping to u, v
	//packed DMs
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &qPack); CHKERRQ(ierr);
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &lambdaPack); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(lambdaPack, pda); CHKERRQ(ierr);

	//x-velocity
	ierr = PetscMalloc(m*sizeof(*lxu), &lxu); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*lyu), &lyu); CHKERRQ(ierr);
	ierr = PetscMemcpy(lxu, lxp, m*sizeof(*lxu)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lyu, lyp, n*sizeof(*lyu)); CHKERRQ(ierr);

	numX = fluid.nx;
	numY = fluid.ny;
	
	// for u, number of nodes in last proc is (numX/m-1) 
	lxu[m-1]--;
	numX = fluid.nx - 1;
	
	// size (lxu * lyu) on each proc only includes  interior domain	
	ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, numX, numY, m, n, 1, 1, lxu, lyu, &uda); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(qPack, uda); CHKERRQ(ierr);

	//y-velocity
	ierr = PetscMalloc(m*sizeof(*lxv), &lxv); CHKERRQ(ierr);
	ierr = PetscMalloc(n*sizeof(*lyv), &lyv); CHKERRQ(ierr);
	ierr = PetscMemcpy(lxv, lxp, m*sizeof(*lxv)); CHKERRQ(ierr);
	ierr = PetscMemcpy(lyv, lyp, n*sizeof(*lyv)); CHKERRQ(ierr);

	numX = fluid.nx;
	numY = fluid.ny;

	//for v, number of nodes in last proc is y direction is (numY/n-1)
	lyv[n-1]--;
	numY = fluid.ny - 1;

	//size (lxv * lyv)  is on each proc, and (numX * numY) is the whole domain(only interior)
	ierr = DMDACreate2d(PETSC_COMM_WORLD, DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_GHOSTED, DMDA_STENCIL_BOX, numX, numY, m, n, 1, 1, lxv, lyv, &vda); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(qPack, vda); CHKERRQ(ierr);

	PetscFree(lxu);
	PetscFree(lyu);
	PetscFree(lxv);
	PetscFree(lyv);

	return 0;
}


 

