#include "NavierStokesSolver.h"
#include <petscdmcomposite.h>  // DMCompositeCreate(MPI_Comm, DM*)
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>  


PetscErrorCode NavierStokesSolver::initialize()
{
	PetscErrorCode ierr;
	ierr = PetscLogStagePush(stageInitialize); CHKERRQ(ierr);
	initialfluid();
	ierr = createDMs(); CHKERRQ(ierr);
	ierr = initializeCommon(); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode NavierStokesSolver::initializeCommon()
{
	PetscErrorCode ierr;

	ierr = createVecs(); CHKERRQ(ierr);

	initializeMeshSpacings(); 
	ierr = initializeFluxes(); CHKERRQ(ierr);
	ierr = initializeLambda(); CHKERRQ(ierr);
	ierr = updateBoundaryGhosts(); CHKERRQ(ierr);

	ierr = createLocalToGlobalMappingsFluxes(); CHKERRQ(ierr);
	ierr = createLocalToGlobalMappingsLambda(); CHKERRQ(ierr);

	ierr = generateDiagonalMatrices(); CHKERRQ(ierr);
	ierr = generateA(); CHKERRQ(ierr);
	ierr = generateBNQ(); CHKERRQ(ierr);
	ierr = generateQTBNQ(); CHKERRQ(ierr);
	ierr = createKSPs(); CHKERRQ(ierr);
	ierr = setNullSpace(); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode NavierStokesSolver::finalize()
{
	PetscErrorCode ierr;

	// DMs
	if(pda != PETSC_NULL) {ierr = DMDestroy(&pda); CHKERRQ(ierr);}
	if(uda != PETSC_NULL) {ierr = DMDestroy(&uda); CHKERRQ(ierr);}
	if(vda != PETSC_NULL) {ierr = DMDestroy(&vda); CHKERRQ(ierr);}
	if(qPack != PETSC_NULL) {ierr = DMDestroy(&qPack); CHKERRQ(ierr);}
	if(lambdaPack != PETSC_NULL) {ierr = DMDestroy(&lambdaPack); CHKERRQ(ierr);}

	//Vecs
	if(q != PETSC_NULL) {ierr = VecDestroy(&q); CHKERRQ(ierr);}
	if(qStar != PETSC_NULL) {ierr = VecDestroy(&qStar); CHKERRQ(ierr);}

	if(qxLocal != PETSC_NULL) {ierr = VecDestroy(&qxLocal); CHKERRQ(ierr);}
	if(qyLocal != PETSC_NULL) {ierr = VecDestroy(&qyLocal); CHKERRQ(ierr);}
	
	if(H != PETSC_NULL) {ierr = VecDestroy(&H); CHKERRQ(ierr);}
	if(rn != PETSC_NULL) {ierr = VecDestroy(&rn); CHKERRQ(ierr);}
	if(bc1 != PETSC_NULL) {ierr = VecDestroy(&bc1); CHKERRQ(ierr);}
	if(rhs1 != PETSC_NULL) {ierr = VecDestroy(&rhs1); CHKERRQ(ierr);}
	if(temp != PETSC_NULL) {ierr = VecDestroy(&temp); CHKERRQ(ierr);}
	if(lambda != PETSC_NULL) {ierr = VecDestroy(&lambda); CHKERRQ(ierr);}
	if(bc2 != PETSC_NULL) {ierr = VecDestroy(&bc2); CHKERRQ(ierr);}
	if(rhs2 != PETSC_NULL) {ierr = VecDestroy(&rhs2); CHKERRQ(ierr);}

	if(uMapping != PETSC_NULL) {ierr = VecDestroy(&uMapping); CHKERRQ(ierr);}
	if(vMapping != PETSC_NULL) {ierr = VecDestroy(&vMapping); CHKERRQ(ierr);}
	if(pMapping != PETSC_NULL) {ierr = VecDestroy(&pMapping); CHKERRQ(ierr);}

	if(MHat != PETSC_NULL) {ierr = VecDestroy(&MHat); CHKERRQ(ierr);}
	if(RInv != PETSC_NULL) {ierr = VecDestroy(&RInv); CHKERRQ(ierr);}
	if(BN != PETSC_NULL) {ierr = VecDestroy(&BN); CHKERRQ(ierr);}

	//Mats
	if(A!=PETSC_NULL) {ierr = MatDestroy(&A); CHKERRQ(ierr);}
	if(QT != PETSC_NULL) {ierr = MatDestroy(&QT); CHKERRQ(ierr);}
	if(BNQ != PETSC_NULL) {ierr = MatDestroy(&BNQ); CHKERRQ(ierr);}
	if(QTBNQ != PETSC_NULL) {ierr = MatDestroy(&QTBNQ); CHKERRQ(ierr);}

	//KSPs
	if(ksp1 != PETSC_NULL) {ierr = KSPDestroy(&ksp1); CHKERRQ(ierr);}
	if(ksp2 != PETSC_NULL) {ierr = KSPDestroy(&ksp2); CHKERRQ(ierr);}

	return 0;
}

PetscErrorCode NavierStokesSolver::generateRHS1()
{
	PetscErrorCode ierr;
	ierr = VecWAXPY(rhs1, 1.0, rn, bc1); CHKERRQ(ierr);
	ierr = VecPointwiseMult(rhs1, MHat, rhs1); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode NavierStokesSolver::generateRHS2()
{
	PetscErrorCode ierr;
	ierr = VecScale(bc2, -1.0); CHKERRQ(ierr);
	ierr = MatMultAdd(QT, qStar, bc2, rhs2); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode NavierStokesSolver::stepTime()
{	
	PetscErrorCode ierr;

	//solve intermediate velocity A q_star = r1
	ierr = PetscLogStagePush(stageSolveIntermediateVelocity); CHKERRQ(ierr);
	ierr = calculateExplicitTerms(); CHKERRQ(ierr);
	ierr = updateBoundaryGhosts(); CHKERRQ(ierr);
	ierr = generateBC1(); CHKERRQ(ierr);
	ierr = generateRHS1();CHKERRQ(ierr);
	ierr = solveIntermediateVelocity();CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	//solve Poisson for pressure
	ierr = PetscLogStagePush(stageSolvePoissonSystem); CHKERRQ(ierr);
	ierr = generateBC2(); CHKERRQ(ierr);
	ierr = generateRHS2(); CHKERRQ(ierr);
	ierr = solvePoissonSystem(); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	//project pressure to satisfy continuity
	
	ierr = PetscLogStagePush(stageProjectionStep); CHKERRQ(ierr);
	ierr = projectionStep(); CHKERRQ(ierr);
	ierr = PetscLogStagePop(); CHKERRQ(ierr);

	timeStep++;

	return 0;
}

PetscErrorCode NavierStokesSolver::solveIntermediateVelocity()
{
	PetscErrorCode ierr; 
	ierr = KSPSolve(ksp1, rhs1, qStar); CHKERRQ(ierr);
	ierr = KSPView(ksp1, PETSC_VIEWER_STDOUT_WORLD);
	return 0;
}

PetscErrorCode NavierStokesSolver::solvePoissonSystem()
{
	PetscErrorCode ierr;
	ierr = KSPSolve(ksp2, rhs2, lambda); CHKERRQ(ierr);
	ierr = KSPView(ksp2,PETSC_VIEWER_STDOUT_WORLD);
	return 0;
}

PetscErrorCode NavierStokesSolver::projectionStep()
{
	// q = q^* - B^N Q \lambda
	PetscErrorCode ierr;
	ierr = MatMult(BNQ, lambda, temp); CHKERRQ(ierr);
	ierr = VecWAXPY(q, -1.0, temp, qStar); CHKERRQ(ierr);
	return 0;
}

/// ?
PetscBool NavierStokesSolver::savePoint()
{
	return (timeStep%fluid.nsave == 0)?PETSC_TRUE:PETSC_FALSE;
}

PetscBool NavierStokesSolver::finished()
{
	return (timeStep >= fluid.nt) ? PETSC_TRUE : PETSC_FALSE;
}

PetscErrorCode NavierStokesSolver:: generateQTBNQ()
{
	PetscErrorCode ierr;
	PetscLogEvent GENERATE_QTBNQ;
	ierr = PetscLogEventRegister("generateQTBNQ",0, &GENERATE_QTBNQ); CHKERRQ(ierr);
	ierr = PetscLogEventBegin(GENERATE_QTBNQ, 0, 0, 0, 0); CHKERRQ(ierr);
	ierr = MatMatMult(QT, BNQ, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &QTBNQ); CHKERRQ(ierr);
	ierr = PetscLogEventEnd(GENERATE_QTBNQ, 0, 0, 0,0); CHKERRQ(ierr);

	return 0;
}

/*
 * cols, array of column indices where non-zeros are present in a particular row of a matrix
 * numCols, number of cols
 * rowStart, start index of the portion of result vector that resident on the current process
 * rowEnd
 * d_nnz number of non_zeros in the diagonal portion of the matrix
 * o_nnz number of non_zeros in the off-diagonal portion of hte matrix
 */

void NavierStokesSolver::countNumNonZeros(PetscInt *cols,size_t numCols, PetscInt rowStart, PetscInt rowEnd, PetscInt &d_nnz, PetscInt &o_nnz)
{
	d_nnz = 0; 
	o_nnz = 0; 
	for(size_t i=0; i<numCols; i++)
	{
		(cols[i] >= rowStart && cols[i] < rowEnd) ? d_nnz++ : o_nnz++;
	}
}

#include "NavierStokes/calculateExplicitTerms.cpp"
#include "NavierStokes/createDMs.cpp"
#include "NavierStokes/createLocalToGlobalMappingsFluxes.cpp"
#include "NavierStokes/createLocalToGlobalMappingsLambda.cpp"
#include "NavierStokes/createVecs.cpp"
#include "NavierStokes/generateA.cpp"
#include "NavierStokes/generateBNQ.cpp"
#include "NavierStokes/generateDiagonalMatrices.cpp"
#include "NavierStokes/initialFluxes.cpp"
#include "NavierStokes/initializeLambda.cpp"
#include "NavierStokes/initializeMeshSpacings.cpp"
#include "NavierStokes/setNullSpace.cpp"
#include "NavierStokes/updateBoundaryGhosts.cpp"
#include "NavierStokes/writeData.cpp"
#include "NavierStokes/writeGrid.cpp"
#include "NavierStokes/createKSPs.cpp"
#include "NavierStokes/generateBC1.cpp"
#include "NavierStokes/generateBC2.cpp"
#include "NavierStokes/initialfluid.cpp"
