#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::writeGrid()
{
	PetscErrorCode	ierr;
	PetscInt	rank;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if(rank==0)
	{
		std::ofstream f("cases/grid.txt");
		for(std::vector<PetscReal>::const_iterator i=fluid.x.begin(); i != fluid.x.end(); i++)
			f << *i << '\n';
		for(std::vector<PetscReal>::const_iterator i=fluid.y.begin(); i != fluid.y.end(); i++)
			f << *i << '\n';
		f.close();
	}

	return 0;
}
