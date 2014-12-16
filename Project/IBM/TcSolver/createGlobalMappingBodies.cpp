#include "../TcSolver.h"

PetscErrorCode	TcSolver::createGlobalMappingBodies()
{
	PetscErrorCode ierr;
	PetscInt	numProcs;
	PetscInt	globalIndex;

	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

	globalIndex = 0;
	for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
	{
		// this globalIndex is used for (pressure and body force)
		globalIndex += numPhiOnProcess[procIdx]; // in each proc, the first numPhi index is for phi, the following index is for body force, and which has two components, +2
		for(std::vector<PetscInt>::iterator i=boundaryPointIndices[procIdx].begin(); i!=boundaryPointIndices[procIdx].end();i++)
		{
			globalIndexMapping[*i] = globalIndex;
			globalIndex += 2;
		}
	}

	return 0;
}

	
