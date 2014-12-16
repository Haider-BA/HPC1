#if !defined(TC_SOLVER_H)
#define TC_SOLVER_H

#include "INS2/NavierStokesSolver.h"

class TcSolver:public NavierStokesSolver
{
	public:
		PetscInt	startGlobalIndex; //
		DM	bda; //body force 
		Mat	ET;
		PetscReal	force[2];
		Vec	nullSpaceVec, regularizedForce;
		//std::ofstream	forceFile;

		std::vector<PetscReal> x,y;
		std::vector<PetscInt>	I, J;
		std::vector<PetscInt>	globalIndexMapping;
		std::vector<PetscInt>	numBoundaryPointsOnProcess;
		std::vector<PetscInt>	numPhiOnProcess;
		std::vector<std::vector<PetscInt> > boundaryPointIndices;

		PetscErrorCode	initializeLambda();
		PetscErrorCode	initializeBodies();
		PetscErrorCode	generateBodyInfo();
		PetscErrorCode	calculateCellIndices();
		PetscErrorCode	createDMs();
		PetscErrorCode	createVecs();
		PetscErrorCode	setNullSpace();
		PetscErrorCode	generateBNQ();
		PetscErrorCode	generateBC2();
		PetscErrorCode	createGlobalMappingBodies();
//		PetscErrorCode	calculateForce();
//		PetscErrorCode	writeForces();
//		PetscErrorCode	writeLambda();

		PetscReal	dhRoma(PetscReal x, PetscReal h);
		PetscReal	delta(PetscReal	x, PetscReal y, PetscReal h);
		PetscBool	isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal xBody, PetscReal yBody, PetscReal radius, PetscReal* delta);

	public:
		PetscErrorCode initialize();
		PetscErrorCode finalize();
		PetscErrorCode	writeData();

		TcSolver()
//son_class inherts bash_class constructor 
	{
		bda = PETSC_NULL;
		ET = PETSC_NULL;
		nullSpaceVec = PETSC_NULL;
		regularizedForce = PETSC_NULL;
	}

};

#endif
