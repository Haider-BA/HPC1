#include "../NavierStokesSolver.h"
#include <iostream>

PetscErrorCode NavierStokesSolver::writeData()
{
	PetscErrorCode ierr;
	PetscInt	rank;
	Vec		qxGlobal, qyGlobal, phi;
	PetscViewer	viewer;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	if(rank==0)
	{
		PetscInt its1, its2;
		ierr = KSPGetIterationNumber(ksp1,&its1);
		ierr = KSPGetIterationNumber(ksp2,&its2);
		std::cout << "Iteration in ksp1" << '\t' << its1 << "Iteration in ksp2" << its2 << std::endl;
	}

     if(timeStep%fluid.nsave ==0)
     {
	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "xflux.dat", FILE_MODE_WRITE, &viewer);
	ierr = VecView(qxGlobal, viewer);
	ierr = PetscViewerDestroy(&viewer);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "yflux.dat", FILE_MODE_WRITE, &viewer);
	ierr = VecView(qyGlobal, viewer);
	ierr = PetscViewerDestroy(&viewer);

	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal);
	
	ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "phi.dat", FILE_MODE_WRITE, &viewer);
	ierr = VecView(phi, viewer);
	ierr = PetscViewerDestroy(&viewer);
	
	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "data written done\n");
     }

     return 0;
}
