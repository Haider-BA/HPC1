#include "../NavierStokesSolver.h"

PetscErrorCode NavierStokesSolver::createVecs()
{
	PetscErrorCode ierr;

	ierr = DMCreateLocalVector(uda, &qxLocal); CHKERRQ(ierr);
	ierr = DMCreateLocalVector(vda, &qyLocal); CHKERRQ(ierr);

	//global vectors
	ierr = DMCreateGlobalVector(qPack, &q); CHKERRQ(ierr); //velocity fluxes
        ierr = VecDuplicate(q, &qStar); CHKERRQ(ierr); //intermediate velocity
	ierr = VecDuplicate(q, &H);	CHKERRQ(ierr); //convective term
	ierr = VecDuplicate(q, &rn);	CHKERRQ(ierr); //explicit term
	ierr = VecDuplicate(q, &bc1);	CHKERRQ(ierr); //B.C. from implicit terms
	ierr = VecDuplicate(q, &rhs1);	CHKERRQ(ierr); //rhs for intermediate solver
	ierr = VecDuplicate(q, &MHat);	CHKERRQ(ierr); // dialogal matrix
	ierr = VecDuplicate(q, &RInv);  CHKERRQ(ierr); // translation form velocity to flux
	ierr = VecDuplicate(q, &BN);	CHKERRQ(ierr); //approxiamte inverse of A
	ierr = VecDuplicate(q, &temp);	CHKERRQ(ierr); 	// ?

	ierr = DMCreateGlobalVector(lambdaPack, &lambda); CHKERRQ(ierr); //pressure
	ierr = VecDuplicate(lambda, &bc2); CHKERRQ(ierr);
	ierr = VecDuplicate(lambda, &rhs2); CHKERRQ(ierr);

	return 0;
}
	
