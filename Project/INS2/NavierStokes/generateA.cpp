#include "../NavierStokesSolver.h"


void getColumns(PetscReal **globalIndices, PetscInt i, PetscInt j, PetscInt *cols)
{
	cols[0] = globalIndices[j][i];
	cols[1] = globalIndices[j][i-1];
	cols[2] = globalIndices[j][i+1];
	cols[3] = globalIndices[j-1][i];	
	cols[4] = globalIndices[j+1][i];
}

void getCoefficients(PetscReal dxMinus, PetscReal dxPlus, PetscReal dyMinus, PetscReal dyPlus, PetscReal *values)
{
	values[0] = -(2.0/dxMinus/dxPlus + 2.0/dyMinus/dyPlus);
	values[1] = 2.0/dxMinus/(dxMinus + dxPlus);
	values[2] = 2.0/ dxPlus/(dxMinus + dxPlus);
	values[3] = 2.0/dyMinus/(dyMinus + dyPlus);
	values[4] = 2.0/ dyPlus/(dyMinus + dyPlus);
}

PetscErrorCode NavierStokesSolver::generateA()
{
	PetscErrorCode ierr;
	PetscInt       mstart, nstart, m, n, i, j;
	PetscInt       cols[5];
	PetscReal      values[5];
	PetscReal      **uGlobalIdx, **vGlobalIdx;
	PetscInt       qStart, qEnd, qLocalSize;
	PetscInt       *d_nnz, *o_nnz;
	PetscInt       localIdx;
	PetscReal      dt = fluid.dt,
	               nu = fluid.nu,
	               alphaImplicit = 0.5;

	// ownership range of q
	ierr = VecGetOwnershipRange(q, &qStart, &qEnd); CHKERRQ(ierr);
	qLocalSize = qEnd-qStart;

	// create arrays to store nnz values
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &d_nnz); CHKERRQ(ierr);
	ierr = PetscMalloc(qLocalSize*sizeof(PetscInt), &o_nnz); CHKERRQ(ierr);

	// determine the number of non-zeros in each row
	// in the diagonal and off-diagonal portions of the matrix
	localIdx = 0;
	// U
	ierr = DMDAVecGetArray(uda, uMapping, &uGlobalIdx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			getColumns(uGlobalIdx, i, j, cols);
			countNumNonZeros(cols, 5, qStart, qEnd, d_nnz[localIdx], o_nnz[localIdx]);
			localIdx++;
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &uGlobalIdx); CHKERRQ(ierr);
	// V
	ierr = DMDAVecGetArray(vda, vMapping, &vGlobalIdx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			getColumns(vGlobalIdx, i, j, cols);
			countNumNonZeros(cols, 5, qStart, qEnd, d_nnz[localIdx], o_nnz[localIdx]);
			localIdx++;
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &vGlobalIdx); CHKERRQ(ierr);

	// create and allocate memory for matrix A
	ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
	ierr = MatSetSizes(A, qLocalSize, qLocalSize, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr);
	ierr = MatSetFromOptions(A); CHKERRQ(ierr);
	ierr = MatSeqAIJSetPreallocation(A, 0, d_nnz); CHKERRQ(ierr);
	ierr = MatMPIAIJSetPreallocation(A, 0, d_nnz, 0, o_nnz); CHKERRQ(ierr);

	// deallocate d_nnz and o_nnz
	ierr = PetscFree(d_nnz); CHKERRQ(ierr);
	ierr = PetscFree(o_nnz); CHKERRQ(ierr);

	// assemble matrix A
	// U
	ierr = DMDAVecGetArray(uda, uMapping, &uGlobalIdx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(uda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			getCoefficients(dxU[i], dxU[i+1], dyU[j], dyU[j+1], values);
			getColumns(uGlobalIdx, i, j, cols);
	
		/*   index i for the No. of rows
		     index j  for the No. of cols
		*/ 
			ierr = MatSetValues(A, 1, &cols[0], 5, cols, values, INSERT_VALUES); CHKERRQ(ierr);	
		}
	}
	ierr = DMDAVecRestoreArray(uda, uMapping, &uGlobalIdx); CHKERRQ(ierr);
	// V
	ierr = DMDAVecGetArray(vda, vMapping, &vGlobalIdx); CHKERRQ(ierr);
	ierr = DMDAGetCorners(vda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);
	for(j=nstart; j<nstart+n; j++)
	{
		for(i=mstart; i<mstart+m; i++)
		{
			getCoefficients(dxV[i], dxV[i+1], dyV[j], dyV[j+1], values);
			getColumns(vGlobalIdx, i, j, cols);
			ierr = MatSetValues(A, 1, &cols[0], 5, cols, values, INSERT_VALUES); CHKERRQ(ierr);
		}
	}
	ierr = DMDAVecRestoreArray(vda, vMapping, &vGlobalIdx); CHKERRQ(ierr);

	ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
	ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

	// for check 
//	ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);

	ierr = MatScale(A, nu*alphaImplicit); CHKERRQ(ierr);
	ierr = MatShift(A, -1.0/dt); CHKERRQ(ierr);
	ierr = MatScale(A, -1.0); CHKERRQ(ierr);
	ierr = MatDiagonalScale(A, MHat, RInv);

	return 0;
}
