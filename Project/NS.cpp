

	
stepTime()
{	
	//solve intermediate velocity A q_star = r1
	petscLogStagePush(satgeSolveIntermediateVelocity);
	calculateExplicitTerms();
	updateBoundaryGhosts();
	generateBC1();
	generateRHS1();
	solveIntermediateVelocity();
	PetsLogStagePop();


	//solve Poisson system for pressure
	PetscLogStagePush;
	generateR2();
	generateRHS2();	
	solvePoissonSystem();
	PetscLogStagePop();

	//project pressure field to satisfy continuity and body force to satisy the no-slip condition
	Push;
	projectionStep();
	Pop;

}

solveIntermediatVelocity()
{
	KSPSolve(ksp1, rhs1, qStar);
}

solvePoissonSystem()
{
	KSPSolve(ksp2, rhs2, lambda);
}

projectionStep()
{
	// q = q^* - B^N Q \lambda
	MatMult(BNQ, lambda, temp);
	VecWAXPY(q, -1.0, temp, qStar);
}

generateQTBNQ()
{
	PetscLogEvent GENERATE_QTBNQ;
	PetscLogEventRegister;
	PetscLogEventBegin();
	MatMatMult(QT, BNQ, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &QTBNQ);
	PetscLogEventEnd();

}

/*
 * cols, array of column indices where non-zeros are present in a particular row of a matrix
 * numCols, number of cols
 * rowStart, start index of the portion of result vector that resident on the current process
 * rowEnd
 * d_nnz number of non_zeros in the diagonal portion of the matrix
 * o_nnz number of non_zeros in the off-diagonal portion of hte matrix
 */

template<PetscInt dim>
void NavierStokesSolver<dim>::countNumNonZeros(*cols, numCols, rowStart, rowEnd, &d_nnz, &o_nnz);



