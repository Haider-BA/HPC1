#include "../TcSolver.h"

PetscErrorCode	TcSolver::generateBNQ()
{
	PetscErrorCode ierr;
	PetscInt	numProcs;
	PetscInt	mstart, nstart, m, n;
	PetscInt	*BNQ_d_nnz, *BNQ_o_nnz;
	PetscInt	*ET_d_nnz, *ET_o_nnz;
	PetscInt	qStart, qEnd, qLocalSize;
	PetscInt	lambdaStart, lambdaEnd, lambdaLocalSize;
	PetscInt	fStart, fEnd, fLocalSize;
	PetscInt	localIdx;
	PetscReal	**pGlobalIdx;
	PetscInt	row, cols[2], BNQ_col, ET_col, value;
	PetscReal	values[2] = {-1.0, 1.0};
	PetscReal	disp[2];
	PetscReal	xCoord, yCoord, h;
	Vec		fGlobal;
	PetscLogEvent	GENERATE_BNQ;

	ierr = PetscLogEventRegister("generateBNQ",0,&GENERATE_BNQ); CHKERRQ(ierr);
	ierr = PetscLogEventBegin(GENERATE_BNQ,0,0,0,0);CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD, &numProcs); CHKERRQ(ierr);

	ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
	qLocalSize = qEnd - qStart;

	//create arrays to store nnz values
	//BNQ
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &BNQ_d_nnz); CHKERRQ(ierr);
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &BNQ_o_nnz); CHKERRQ(ierr);
	//ET
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &ET_d_nnz);CHKERRQ(ierr);
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt),&ET_o_nnz); CHKERRQ(ierr);

	ierr = VecGetOwnershipRange(lambda, &lambdaStart, &lambdaEnd); CHKERRQ(ierr);
	lambdaLocalSize =  lambdaEnd - lambdaStart;

	//ownership range of f
	ierr = DMCompositeGetAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);
	ierr = VecGetOwnershipRange(fGlobal, &fStart, &fEnd); CHKERRQ(ierr);
	fLocalSize = fEnd - fStart;
	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, NULL, &fGlobal); CHKERRQ(ierr);

	//get mapping of pressure values
	ierr = DMDAVecGetArray(pda, pMapping, &pGlobalIdx); CHKERRQ(ierr);

	//update position 
//	ierr = updateBody();

	//determine the number of non-zeros in each row
	//in the diagonal and off diagonal portions of matrix
	localIdx = 0;
	//U 
	ierr = DMDAGetCorners(uda,&mstart,&nstart,NULL,&m,&n,NULL); CHKERRQ(ierr);
	for(PetscInt j=nstart; j<nstart+n;j++)
	{
		yCoord = 0.5*(fluid.y[j] + fluid.y[j+1]);
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			h=fluid.dx;
			xCoord=fluid.x[i+1];
			//G portion
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j][i+1];
			countNumNonZeros(cols,2,lambdaStart,lambdaEnd,BNQ_d_nnz[localIdx],BNQ_o_nnz[localIdx]);
			//ET portion
			ET_d_nnz[localIdx] = 0;
			ET_o_nnz[localIdx] = 0;
			PetscInt numPhi = 0;
			for(PetscInt procIdx = 0; procIdx<numProcs; procIdx++)
			{
				numPhi += numPhiOnProcess[procIdx];
				for( std::vector<PetscInt>::iterator l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end(); l++)
				{
					//current cell and its neigbhours searching
					if(j >= I[*l]-2 && i<= I[*l]+2 && j >= J[*l]-2 && j<= J[*l]+2)
					{
						if(isInfluenced(xCoord, yCoord, x[*l], y[*l], 1.5*h, disp))
						{
							BNQ_col = globalIndexMapping[*l];
							(BNQ_col >= lambdaStart && BNQ_col < lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
							ET_col = globalIndexMapping[*l] - numPhi;
							(ET_col >= fStart && ET_col <fEnd)? ET_d_nnz[localIdx]++ : ET_o_nnz[localIdx]++;
						}
					}
				}
			}
			localIdx++;
		}
	}

	//V
	ierr = DMDAGetCorners(vda,&mstart,&nstart,NULL,&m,&n,NULL); CHKERRQ(ierr);
	for(PetscInt j=nstart; j<nstart+n;j++)
	{
		h = fluid.dy;
		yCoord = fluid.y[j+1];
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			xCoord=0.5*(fluid.x[i] + fluid.x[i+1]);
			//G portion
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j+1][i];
			countNumNonZeros(cols,2,lambdaStart,lambdaEnd,BNQ_d_nnz[localIdx], BNQ_o_nnz[localIdx]);
			//ET portion
			ET_d_nnz[localIdx] = 0;
			ET_o_nnz[localIdx] = 0;
			PetscInt numPhi = 0;
			for(PetscInt procIdx=0; procIdx<numProcs; procIdx++)
			{
				numPhi += numPhiOnProcess[procIdx];
				for(std::vector<PetscInt>::iterator l=boundaryPointIndices[procIdx].begin();l!=boundaryPointIndices[procIdx].end(); l++)
				{
					if(i>=I[*l]-2 && i<=I[*l]+2 && j>=J[*l]-2 && j<=J[*l]+2)
					{
						if(isInfluenced(xCoord,yCoord,x[*l],y[*l],1.5*h, disp))
						{
							BNQ_col = globalIndexMapping[*l]+1;
							(BNQ_col>=lambdaStart && BNQ_col <lambdaEnd) ? BNQ_d_nnz[localIdx]++ : BNQ_o_nnz[localIdx]++;
							ET_col = globalIndexMapping[*l] - numPhi + 1;
							(ET_col >= fStart && ET_col < fEnd) ? ET_d_nnz[localIdx]++:ET_o_nnz[localIdx]++;
						}
					}
				}
			}
			localIdx++;
		}
	}

	//allocate memory for the matrices
	//BNQ
	ierr = MatCreate(PETSC_COMM_WORLD,&BNQ);CHKERRQ(ierr);
	ierr = MatSetSizes(BNQ,qLocalSize,lambdaLocalSize,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
	ierr = MatSetFromOptions(BNQ); CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(BNQ,0,BNQ_d_nnz);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(BNQ,0,BNQ_d_nnz,0,BNQ_o_nnz);CHKERRQ(ierr);
	//ET
	ierr = MatCreate(PETSC_COMM_WORLD,&ET);CHKERRQ(ierr);
	ierr = MatSetSizes(ET, qLocalSize, fLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
	ierr = MatSetFromOptions(ET); CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(ET,0,ET_d_nnz);CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(ET,0,ET_d_nnz,0,ET_o_nnz);CHKERRQ(ierr);

	//deallocate nnz arrays
	//BNQ
	ierr = PetscFree(BNQ_d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(BNQ_o_nnz); CHKERRQ(ierr);
	//ET
	ierr = PetscFree(ET_d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(ET_o_nnz); CHKERRQ(ierr);

	//assemble matrices Q and ET
	localIdx = 0;
	//U
	ierr = DMDAGetCorners(uda,&mstart,&nstart,NULL,&m,&n,NULL);CHKERRQ(ierr);
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		yCoord = 0.5*(fluid.y[j]+fluid.y[j+1]);
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			h = fluid.dx;
			xCoord = fluid.x[i+1];
			row = localIdx + qStart;
			//G portion
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j][i+1];
			ierr = MatSetValues(BNQ,1,&row,2,cols,values,INSERT_VALUES);CHKERRQ(ierr);
			//ET portion
			PetscInt	numPhi = 0;
			for(PetscInt procIdx = 0; procIdx<numProcs; procIdx++)
			{
				numPhi += numPhiOnProcess[procIdx];
				for(std::vector<PetscInt>::iterator l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end();l++)
				{
					if(i>=I[*l]-2 && i<=I[*l]+2 && j>=J[*l]-2 && j<=J[*l]+2)
					{
						if(isInfluenced(xCoord,yCoord,x[*l], y[*l],1.5*h, disp))
						{
							BNQ_col = globalIndexMapping[*l];
							value = h*delta(disp[0],disp[1],h);
							ierr = MatSetValue(BNQ,row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
							ET_col = globalIndexMapping[*l] - numPhi;
							ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
						}
					}
				}
			}
		localIdx++;
		}
	}
	//V
	ierr = DMDAGetCorners(vda,&mstart,&nstart,NULL,&m,&n,NULL);CHKERRQ(ierr);
	for(PetscInt j=nstart; j<nstart+n; j++)
	{
		h = fluid.y[j];
		yCoord = fluid.y[j+1];
		for(PetscInt i=mstart; i<mstart+m; i++)
		{
			xCoord = 0.5*(fluid.x[i] +  fluid.x[i+1]);
			row = localIdx + qStart;
			//G portion
			cols[0] = pGlobalIdx[j][i];
			cols[1] = pGlobalIdx[j+1][i];
			ierr = MatSetValues(BNQ,1,&row,2,cols,values,INSERT_VALUES);CHKERRQ(ierr);
			//ET portion
			PetscInt	numPhi = 0;
			for(PetscInt procIdx = 0; procIdx<numProcs; procIdx++)
			{
				numPhi += numPhiOnProcess[procIdx];
				for(std::vector<PetscInt>::iterator l=boundaryPointIndices[procIdx].begin(); l!=boundaryPointIndices[procIdx].end();l++)
				{
					if(i>=I[*l]-2 && i<=I[*l]+2 && j>=J[*l]-2 && j<=J[*l]+2)
					{
						if(isInfluenced(xCoord,yCoord,x[*l], y[*l],1.5*h, disp))
						{
							BNQ_col = globalIndexMapping[*l];
							value = h*delta(disp[0],disp[1],h);
							ierr = MatSetValue(BNQ,row, BNQ_col, value, INSERT_VALUES); CHKERRQ(ierr);
							ET_col = globalIndexMapping[*l] - numPhi;
							ierr = MatSetValue(ET, row, ET_col, value, INSERT_VALUES); CHKERRQ(ierr);
						}
					}
				}
			}
		localIdx++;
		}
	}

		ierr = DMDAVecRestoreArray(pda,pMapping,&pGlobalIdx);CHKERRQ(ierr);
		//assemble the matrices
		//BNQ
		ierr = MatAssemblyBegin(BNQ,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
		ierr = MatAssemblyEnd(BNQ,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
		//ET
		ierr = MatAssemblyBegin(ET,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(ET,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

		ierr = MatTranspose(BNQ,MAT_INITIAL_MATRIX,&QT); CHKERRQ(ierr);
		ierr = MatDiagonalScale(BNQ,BN,NULL); CHKERRQ(ierr);

		ierr = PetscLogEventEnd(GENERATE_BNQ,0,0,0,0); CHKERRQ(ierr);

		return 0;
}


