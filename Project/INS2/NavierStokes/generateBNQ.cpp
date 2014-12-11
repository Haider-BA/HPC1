#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::generateBNQ()
{
	PetscErrorCode ierr;
	PetscInt	i,j;
	PetscInt	mstart, nstart, m,n;
	PetscInt	*d_nnz, *o_nnz;
	PetscInt	qStart, qEnd, lambdaStart, lambdaEnd, qLocalSize, lambdaLocalSize;
	PetscInt	localIdx;
	PetscReal	**pGlobalIdx;
	PetscInt	row, cols[2];
	PetscReal	values[2] = {-1.0, 1.0};

	ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
	qLocalSize = qEnd - qStart;

	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); 
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz);

	ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
	lambdaLocalSize = lambdaEnd - lambdaStart;

	ierr = DMDAVecGetArray(pda,pMapping,&pGlobalIdx);
	
	//determine number of non-zeros in each row in the diagonal and off-diagonal portions of the matrix
	localIdx = 0;
	//U
	ierr = DMDAGetCorners(uda,&mstart, &nstart, NULL, &m, &n, NULL); 
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart;i<mstart+m;i++)
		{
			cols[0] = pGlobalIdx[j][i]; //current position
			cols[1] = pGlobalIdx[j][i+1]; //current right position(along x+ direction)
			countNumNonZeros(cols,2,lambdaStart,lambdaEnd,d_nnz[localIdx], o_nnz[localIdx]);
			localIdx++;
		}
	}
	//V
	ierr = DMDAGetCorners(vda,&mstart,&nstart,NULL,&m,&n,NULL);
	for(j=nstart;j<nstart+n;j++)
	{
		for(i=mstart;i<mstart+m;i++)
		{
			cols[0] = pGlobalIdx[j][i]; 	
			cols[1] = pGlobalIdx[j+1][i]; //current top position(along y+ direction)
			countNumNonZeros(cols,2,lambdaStart,lambdaEnd,d_nnz[localIdx],o_nnz[localIdx]);
			localIdx++;
		}
	}

	//allocate memory for the matrix
	ierr = MatCreate(PETSC_COMM_WORLD, &BNQ); CHKERRQ(ierr);
	ierr = MatSetSizes(BNQ, qLocalSize, lambdaLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
	ierr = MatSetFromOptions(BNQ);
	ierr = MatSeqAIJSetPreallocation(BNQ,0,d_nnz);
	ierr = MatMPIAIJSetPreallocation(BNQ,0,d_nnz,0,o_nnz);

	ierr = PetscFree(d_nnz);
	ierr = PetscFree(o_nnz);

	//assemble matrix Q
	localIdx = 0;
	//U
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m;i++)
		{
			row = localIdx + qStart;
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j][i+1];
			ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
			localIdx++;
		}
	}

	//V
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m;i++)
		{
			row = localIdx + qStart;
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j][i+1];
			ierr = MatSetValues(BNQ, 1, &row, 2, cols, values, INSERT_VALUES); CHKERRQ(ierr);
			localIdx++;
		}
	}

	ierr = DMDAVecRestoreArray(pda, pMapping, &pGlobalIdx); CHKERRQ(ierr);
	ierr = MatAssemblyBegin(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(BNQ, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	//for check
//	PetscPrintf(MPI_COMM_WORLD,"Q matrix\n");
//	MatView(BNQ, PETSC_VIEWER_STDOUT_WORLD);

	ierr = MatTranspose(BNQ, MAT_INITIAL_MATRIX, &QT); CHKERRQ(ierr);
	ierr = MatDiagonalScale(BNQ, BN, NULL); CHKERRQ(ierr);

//	PetscPrintf(MPI_COMM_WORLD,"BNQ matrix\n");
//	MatView(BNQ, PETSC_VIEWER_STDOUT_WORLD);

	return 0;
}
