#if !defined(NAVIER_STOKES_SOLVER_H)
#define	NAVIER_STOKES_SOLVER_H

#include <petscdmda.h>
#include <petscksp.h>
#include <fstream>
#include <vector>


struct fluidsimulation
 {
		PetscReal nu; //kinematic viscosity 
		PetscReal initialVelocity[2];
		PetscReal bc[2][4]; // u, v-flux at four walls

		PetscReal	dt; //timestep
		PetscInt	nt, //number of timestep
				nsave; //intervals at which simulaiton data is saved
		
		PetscInt	nx, //number of cells in x-direction
				ny;

		PetscReal	dx, dy; //cell-spacing of x, y
		std::vector<PetscReal> x, //x-coordinates of nodes
				       y; 
	
};


class NavierStokesSolver
{
public:
	fluidsimulation  fluid; 

	PetscInt	timeStep,
			iterationCount1,
			iterationCount2;

	std::vector<PetscReal> dxU, dyU, 
			       dxV, dyV;

	DM	pda,
		uda,
		vda,
		qPack,
		lambdaPack;

	Vec	qxLocal,
		qyLocal;

	Vec	uMapping,
		vMapping,
		pMapping;

	Vec	H, rn;
	Vec	RInv, MHat;

	Mat	A;
	Mat	QT, BNQ;
	Mat	QTBNQ;
	Vec	BN;
	Vec	bc1, rhs1, bc2, rhs2, temp;
	Vec	q, qStar, lambda;
	KSP	ksp1, ksp2;
	PC	pc2;

	PetscLogStage	stageInitialize,
			stageSolveIntermediateVelocity,
			stageSolvePoissonSystem,
			stageProjectionStep;

PetscErrorCode initializeCommon(); // common data initialize

PetscErrorCode	createDMs(); //create DMDA structures for flow variables

PetscErrorCode	createVecs();  //create Vectors to store flow variables

PetscErrorCode	createKSPs();  //set up Kyrlov solvers used to solve the linear systems

void initializeMeshSpacings();  //initialize the spaces between adjacent velocity nodes

void  initialfluid(); //initial fluid paramters and simulation paramters

PetscErrorCode initializeFluxes(); //populate flux vector with initial conditions
		
PetscErrorCode	readFluxes(); //red fluxes from previous saved data

PetscErrorCode	initializeLambda();

PetscErrorCode	createLocalToGlobalMappingsFluxes(); //mapping from local flux variables to global flux vector

PetscErrorCode createLocalToGlobalMappingsLambda();

PetscErrorCode	updateBoundaryGhosts(); // update ghost nodes on the domain boundary

PetscErrorCode	generateDiagonalMatrices(); //generate diagonal matrices (M, Rinv)

/*count number of non_zeros in diagonal and off-diagonal protions of the parallel matrices */
void countNumNonZeros(PetscInt *cols, size_t numCols, PetscInt rawStart, PetscInt rowEnd, PetscInt &d_nnz, PetscInt &o_nnz);

PetscErrorCode generateA(); //Generate matrix A

PetscErrorCode calculateExplicitTerms(); //calculate explicit convection and diffusion terms

/*assemble the vector arising from the B.C. in RHS of intermediate velocity solve  */
PetscErrorCode generateBC1();	

PetscErrorCode generateRHS1(); //calculate RHS of intermediate velocity 

PetscErrorCode generateBC2(); //assemble vector arising from b.c. in RHS of pressure-force solve

PetscErrorCode generateRHS2();

PetscErrorCode generateBNQ(); //assemble matrix BNQ ?? JCP 2007  (26)

PetscErrorCode generateQTBNQ(); //assemble QTBNQ (27)

PetscErrorCode setNullSpace(); //specify to Krylov solver the  null space of the LHS matrix in pressure-force solve

PetscErrorCode solveIntermediateVelocity(); //sovle for intermediate velocity flux

PetscErrorCode solvePoissonSystem(); //solve for pressure and body forces

PetscErrorCode projectionStep(); //project pressure nad forces on ot the velocity field to obtain the velocity at next timestep


public:
	PetscErrorCode initialize(); //initial set-up system
	PetscErrorCode finalize(); 
	PetscErrorCode stepTime(); //move simulation forward on timestep
	PetscErrorCode writeData();
//	PetscErrorCode writeGrid();
	PetscBool savePoint(); 
	PetscBool finished();

	NavierStokesSolver()
   {
	//DMs
	pda = PETSC_NULL;
	uda = PETSC_NULL;
	vda = PETSC_NULL;
	qPack = PETSC_NULL; // JCP 2007 (27)
	lambdaPack = PETSC_NULL; // JCP 2007 (27)

	//Vecs
	qxLocal = PETSC_NULL;
	qyLocal = PETSC_NULL;
	q = PETSC_NULL;
	qStar	= PETSC_NULL;
	H	= PETSC_NULL;
	rn	= PETSC_NULL;
	bc1	= PETSC_NULL;
	rhs1	= PETSC_NULL;
	bc2 	= PETSC_NULL;
	rhs2	= PETSC_NULL;
	temp	= PETSC_NULL;
	RInv	= PETSC_NULL;
	MHat	= PETSC_NULL;
	BN	= PETSC_NULL;
	pMapping = PETSC_NULL;
	uMapping = PETSC_NULL;
	vMapping = PETSC_NULL;

	//mats
	A 	= PETSC_NULL;
	QT	= PETSC_NULL;
	BNQ	= PETSC_NULL;
	QTBNQ	= PETSC_NULL;

	//KSP
	ksp1	= PETSC_NULL;
	ksp2	= PETSC_NULL;

	//PCs
	pc2	= PETSC_NULL;

	//PetscLogStage
	PetscLogStageRegister("initialize", &stageInitialize);
	PetscLogStageRegister("solveIntVel", &stageSolveIntermediateVelocity);
	PetscLogStageRegister("solvePoissSys", &stageSolvePoissonSystem);
	PetscLogStageRegister("projectionStep", &stageProjectionStep);

   }

};

#endif


