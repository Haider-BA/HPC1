#include "../TcSolver.h"

PetscErrorCode	TcSolver::generateBodyInfo()
{
	PetscErrorCode	ierr;
	PetscInt	m,n;
	const	PetscInt *lxp, *lyp;
	PetscInt	numProcs;

	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

/* set arrays to store boundaryPoint by number of procs */
	boundaryPointIndices.resize(numProcs);
	numBoundaryPointsOnProcess.resize(numProcs);
	numPhiOnProcess.resize(numProcs);

	globalIndexMapping.resize(x.size());
	ierr = DMDAGetOwnershipRanges(pda,&lxp,&lyp,NULL);CHKERRQ(ierr);
	ierr = DMDAGetInfo(pda, NULL, NULL, NULL, NULL, &m, &n, NULL, NULL,NULL,NULL, NULL, NULL, NULL); CHKERRQ(ierr);

	PetscInt xStart, yStart, xEnd, yEnd,
		 procIdx = 0;

	yStart=0;
	for(PetscInt j=0; j<n; j++)
	{
		yEnd = yStart + lyp[j];
		xStart = 0;
		for(PetscInt i=0; i<m; i++)
		{
			procIdx = j*m + i;
			xEnd = xStart + lxp[i];
			numPhiOnProcess[procIdx] = lxp[i] * lyp[j];
			for(size_t l=0; l<x.size(); l++)
			{
				if(x[l] >= fluid.x[xStart] && x[l] < fluid.x[xEnd] && y[l] >= fluid.y[yStart] && y[l] < fluid.y[yEnd])
				{
					numBoundaryPointsOnProcess[procIdx]++;
				}
			}

			boundaryPointIndices[procIdx].reserve(numBoundaryPointsOnProcess[procIdx]);
			for(size_t l=0; l<x.size(); l++)
			{
				if(x[l]>=fluid.x[xStart] && x[l]<fluid.x[xEnd] && y[l]>= fluid.y[yStart] && y[l]< fluid.y[yEnd])
				{
			/* */
					boundaryPointIndices[procIdx].push_back(l);
				}
			}
			xStart = xEnd;
		}
		yStart = yEnd;
	}
	return 0;
}

				
