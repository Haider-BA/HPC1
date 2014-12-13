#include "NavierStokesSolver.h"

int main(int argc, char** argv)
{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc, &argv, (char*)0, NULL);
	NavierStokesSolver* solver = new NavierStokesSolver(); 
	ierr = solver->initialize();
	solver->writeGrid();
	
	while(!solver->finished())
	{
		ierr = solver->stepTime();
		ierr = solver->writeData();
	}

	ierr = solver->finalize();

	ierr = PetscFinalize();

	return 0;
}

