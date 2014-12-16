#include "../TcSolver.h"

PetscErrorCode TcSolver::updateBody()
{
	PetscErrorCode ierr;
	PetscInt	rank;
	PetscInt	totalPoints;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if(rank ==0)
	{
		totalPoints = x.size();
		//update body position
		for(PetscInt i=0; i<totalPoints; i++)
		{
			x[i] += velx * simParams->dt;
			y[i] += vely * simParams->dt;
		}
	}

	ierr = MPI_Bcast(&x.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Bcast(&y.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

	return 0;
}
