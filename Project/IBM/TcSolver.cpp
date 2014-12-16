
#include "TcSolver.h"
#include <petscdmcomposite.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/stat.h>

PetscErrorCode TcSolver::initialize()
{
	PetscErrorCode ierr;
	ierr = PetscLogStagePush(NavierStokesSolver::stageInitialize); CHKERRQ(ierr);
	NavierStokesSolver::initialfluid();
	ierr = initializeBodies();CHKERRQ(ierr);
	ierr = calculateCellIndices(); CHKERRQ(ierr);
	ierr = createDMs(); CHKERRQ(ierr);
	ierr = createGlobalMappingBodies(); CHKERRQ(ierr);
	ierr = NavierStokesSolver::initializeCommon();CHKERRQ(ierr); //inhert base_class function
	ierr = PetscLogStagePop();CHKERRQ(ierr);

	return 0;
}

PetscErrorCode TcSolver::finalize()
{
	PetscErrorCode ierr;
	ierr = NavierStokesSolver::finalize();

	if(bda != PETSC_NULL) {ierr = DMDestroy(&bda); CHKERRQ(ierr);}
	
	if(ET != PETSC_NULL) {ierr = MatDestroy(&ET); CHKERRQ(ierr);}

	if(regularizedForce != PETSC_NULL) {ierr = VecDestroy(&regularizedForce); CHKERRQ(ierr);}
	if(nullSpaceVec != PETSC_NULL) {ierr = VecDestroy(&nullSpaceVec);CHKERRQ(ierr);}

	return 0;
}

PetscErrorCode TcSolver::createDMs()
{
	PetscErrorCode ierr;
	ierr = NavierStokesSolver::createDMs(); CHKERRQ(ierr);
	ierr = generateBodyInfo(); CHKERRQ(ierr);
	ierr = DMDACreate1d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE,x.size(), 2, 0, &numBoundaryPointsOnProcess.front(), &bda); CHKERRQ(ierr);
// x.size = number of Lagrangian points
// 2 = number of DoF per node(fx, fy)
	ierr = DMCompositeAddDM(NavierStokesSolver::lambdaPack, bda); CHKERRQ(ierr);

	return 0;
}

PetscErrorCode TcSolver::createVecs()
{
	PetscErrorCode ierr;
	ierr = NavierStokesSolver::createVecs();
	ierr =VecDuplicate(NavierStokesSolver::q, &regularizedForce); CHKERRQ(ierr); //
	ierr = VecDuplicate(NavierStokesSolver::lambda, &nullSpaceVec); CHKERRQ(ierr);

	return 0;
}

/* PetscErrorCode TcSolver::writeData()
{
	NavierStokesSolver::writeData();
	writeData();
	calculateForce();
	writeForces();
	return 0;
}
*/

PetscReal TcSolver::dhRoma(PetscReal x, PetscReal h)
{
	 PetscReal r= fabs(x)/h;
	 if(r>1.5) return 0.0;
	 if(r>0.5 && r<= 1.5) return 1.0/(6*h)*(5.0-3.0*r-sqrt(-3.0*(1-r)*(1-r)+1.0));
	 return 1.0/(3*h)*(1.0+sqrt(-3.0*r*r+1.0));
}

PetscReal TcSolver::delta(PetscReal x, PetscReal y, PetscReal h)
{
	return dhRoma(x,h)*dhRoma(y,h);
}

#include "TcSolver/setNullSpace.cpp"
#include "TcSolver/calculateCellIndices.cpp"
#include "TcSolver/initializeLambda.cpp"
#include "TcSolver/generateBodyInfo.cpp"
#include "TcSolver/generateBNQ.cpp"
#include "TcSolver/generateBC2.cpp"
#include "TcSolver/initializeBodies.cpp"
#include "TcSolver/createGlobalMappingBodies.cpp"
#include "TcSolver/isInfluenced.cpp"
//#include "TcSolver/writeLambda.cpp"
//#include "TcSolver/calculateForce.cpp"
//#include "TcSolver/writeForces.cpp"
#include "TcSolver/writeData.cpp"
