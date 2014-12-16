#include "../TcSolver.h"
PetscErrorCode TcSolver::initializeBodies()
{
	PetscErrorCode ierr;
	PetscInt	rank;
	PetscInt	totalPoints;
	PetscReal	velx, vely;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
	if(rank==0)
	{
		PetscReal cx, cy, R;
		PetscInt	numPoints;
		
		cx = 0.1/2;
		cy = 0.1/2;
		R =  0.02;
		numPoints = 8;
		velx = 0.0;
		vely = 0.01;
		
		x.reserve(numPoints);
		y.reserve(numPoints);
		for(PetscInt i=0; i<numPoints; i++)
		{
			x.push_back(cx + R*cos(i*2*PETSC_PI/numPoints));
			y.push_back(cy + R*sin(i*2*PETSC_PI/numPoints));
		}
		totalPoints = x.size();
	}

	ierr = MPI_Bcast(&totalPoints, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Bcast(&velx, 1, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Bcast(&vely, 1, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

	x.resize(totalPoints);
	y.resize(totalPoints);

	ierr = MPI_Bcast(&x.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
	ierr = MPI_Bcast(&y.front(), totalPoints, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

	return 0;
}
