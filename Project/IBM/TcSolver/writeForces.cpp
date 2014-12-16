#include "../TcSolver.h"
#include <fstream>

PetscErrorCode	TcSolver::writeForces()
{
	PetscErrorCode	ierr;
	PetscInt	rank;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

	if(rank==0)
	{
		std::string filename =  "../cases/forces.txt";
		std::ofstream	forceFile;
		if(timeStep ==1)
		{
			forceFile.open(filename.c_str());
		}
		else
		{
			forceFile.open(filename.c_str(), std::ios::out | std::ios::app);
		}

		forceFile << timeStep*fluid.dt << '\t' << force[0] << '\t' << force[1] << std::endl;
		forceFile.close();
	}
		return 0;
}

