#include "TcSolver.h"

int main(int argc, char** argv)
{
	PetscErrorCode ierr;

	ierr = PetscInitialize(&argc, &argv, (char*)0, NULL);
	TcSolver* solver = new TcSolver(); 
	ierr = solver->initialize();
	solver->writeGrid();
	
	while(!solver->finished())
	{
		ierr = solver->stepTime();
//		ierr = solver->calculateForce();
		ierr = solver->writeData();
	}

	ierr = solver->finalize();

	ierr = PetscFinalize();

	return 0;
}

