#include "../NavierStokesSolver.h"
#include <iostream>
#include <string>
#include <sstream>

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
		std::cout << "Iteration in ksp1" << its1 << "\n Iteration in ksp2" << its2 << "\n"<< std::endl;
	}

     if(timeStep%fluid.nsave ==0)
     {
	std::stringstream file_name;
	if(timeStep<10)
		file_name = "00000" + to_string(timeStep);
	else if(timeStep<100)
		file_name = "0000" + to_string(timeStep);
	else if(timeStep<1000)
		file_name = "000" + to_string(timeStep);
	else
		file_name = to_string(timeStep) + ".dat"; 	

	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,  (file_name+"qx.dat").c_str(), FILE_MODE_WRITE, &viewer);
	ierr = VecView(qxGlobal, viewer);
	ierr = PetscViewerDestroy(&viewer);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, (file_name+"qy.dat").c_str(), FILE_MODE_WRITE, &viewer);
	ierr = VecView(qyGlobal, viewer);
	ierr = PetscViewerDestroy(&viewer);

	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal);
	
	ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, (file_name+"phi.dat").c_str(), FILE_MODE_WRITE, &viewer);
	ierr = VecView(phi, viewer);
	ierr = PetscViewerDestroy(&viewer);
	
	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "data written done\n");
     }

     return 0;
}
