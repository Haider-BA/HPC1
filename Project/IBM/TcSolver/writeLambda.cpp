#include "../TcSolver.h"

PetscErrorCode TcSolver::writeLambda()
{
	PetscErrorCode	ierr;
	Vec	phi, fTilde;
	std::string	savePointDir, fileName;
	PetscViewer	viewer;

	//create output folder
	std::stringstream ss;
	ss <<  "cases/ibm/" << std::setfill('0') << std::setw(7) << NavierStokesSolver::timeStep;
	savePointDir = ss.str();

	ierr = DMCompositeGetAccess(NavierStokesSolver::lambdaPack, NavierStokesSolver::lambda, &phi, &fTilde); CHKERRQ(ierr);

	//print phi to file
	ss.str(" ");
	ss.clear();

	ss << savePointDir << "/phi.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
	ierr = VecView(phi, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr); 

	//print fTilde to file
	ss.str(" ");
	ss.clear();

	ss << savePointDir << "/fTilde.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
	ierr = VecView(fTilde, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr); 

	ierr = DMCompositeRestoreAccess(NavierStokesSolver::lambdaPack, NavierStokesSolver::lambda, &phi, &fTilde); CHKERRQ(ierr);

	return 0;
}


