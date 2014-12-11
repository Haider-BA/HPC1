#include <fstream>

void  NavierStokesSolver::initialfluid()
{
	PetscInt rank, i, j;

	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if(rank==0) //read input file only on process 0
	{
		fluidsimulation temp = {0.025, {1.0, 0.0}, {{1.0, 1.0, 1.0, 1.0}, {0.0, 0.0,0.0,0.0}}, 
						0.01, 5, 5, 10, 10, 0.05, 0.05};
		fluid = temp; 		
	}
		MPI_Barrier(PETSC_COMM_WORLD);

		//broadcast number of cells to all process
		MPI_Bcast(&fluid.nu, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
		MPI_Bcast(fluid.initialVelocity, 2, MPIU_REAL, 0, PETSC_COMM_WORLD);
		MPI_Bcast(fluid.bc[0], 4, MPIU_REAL, 0, PETSC_COMM_WORLD);
		MPI_Bcast(fluid.bc[1], 4, MPIU_REAL, 0, PETSC_COMM_WORLD);

		MPI_Bcast(&fluid.dt, 1,MPIU_REAL, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&fluid.nt, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&fluid.nsave, 1,MPIU_INT, 0,PETSC_COMM_WORLD);
		MPI_Bcast(&fluid.nx, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&fluid.ny, 1, MPIU_INT, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&fluid.dx, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&fluid.dy, 1, MPIU_REAL, 0, PETSC_COMM_WORLD);

		//mesh gridding
		fluid.x.resize(fluid.nx +1);
		fluid.y.resize(fluid.ny +1);

		if(rank==0)
		{
			for( i = 0; i< fluid.nx +1; i++)
			{
				fluid.x[i] = fluid.dx * i;
			}
			for(j=0; j<fluid.ny +1; j++)
			{
				fluid.y[j] = fluid.dy * j;
			}

		}
		MPI_Barrier(PETSC_COMM_WORLD);

		MPI_Bcast(&fluid.x.front(), fluid.nx+1, MPIU_REAL, 0, PETSC_COMM_WORLD);
		MPI_Bcast(&fluid.y.front(), fluid.ny+1, MPIU_REAL, 0, PETSC_COMM_WORLD);

}


