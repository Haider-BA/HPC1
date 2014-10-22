

initializeCommon(); // common data initialize

createDMs(); //create DMDA structures for flow variables

createVecs();  //create Vectors to store flow variables

createKSPs();  //set up Kyrlov solvers used to solve the linear systems

initializedMeshSpacings();  //initialize the spaces between adjacent velocity nodes

initializeFluxes(); //flux vector with initial conditions

readFluxes(); //red fluxes from previous saved data


createLocalToGlobalMappingFluxes(); //mapping from local flux variables to global flux vector

updateBoundaryGhosts(); // update ghost nodes on the domain boundary

generateDiagonalMatrices(); //generate diagonal matrices (M, Rinv)

/*count number of non_zeros in diagonal and off-diagonal protions of the parallel matrices */
countNumNonZeros(*cols, size_t, rawStart, rowEnd, &d_nnz, &o_nnz);

generateA(); //Generate matrix A

calculateExplicitTerms(); //calculate explicit convection and diffusion terms

/*assemble the vector arising from the B.C. in RHS of intermediate velocity solve  */
generateBC1();	

generateRHS1(); //calculate RHS of intermediate velocity 

generateR2(); //assemble vector arising from b.c. in RHS of pressure-force solve

generateRHS2();

generateBNQ(); //assemble matrix BNQ ?? JCP 2007  (26)

generateQTBNQ(); //assemble QTBNQ (27)

virtual setNullSpace(); //specify to Krylov solver the  null space of the LHS matrix in pressure-force solve

solveIntermediateVelocity(); //sovle for intermediate velocity flux

solvePoissonSystem(); //solve for pressure and body forces

projectionStep(); //project pressure nad forces on ot the velocity field to obtain the velocity at next timestep

writeFluxes();

public:
	initialize(); //initial set-up system
	finalize(); 
	stepTime(); //move simulation forward on timestep
	writeData();
	writeSimulationInfo();
	writeGrid();
	savePoint(); 
	finished();

	NavierStokesSolver(string, FD, SP, CM)
{
	//class
	caseFolder
	flowDesc
	simParams
	mesh
	timestep

	//DMs
	pda
	uda vda wda
	qPack; // JCP 2007 (27)
	lambdaPack; // JCP 2007 (27)

	//Vecs
	qxLocal qyLocal qzLocal
	q 
	...

	//mats
	A 
	QT
	BNQ
	QTBNQ

	//KSP
	ksp1
	ksp2

	//PCs
	pc2

	//PetscLogStage
	PetscLogStageRegister("initialize");
PPPetscLogStageRegister("solveIntVel");

PetscLogStageRegister("solvePoissSys");

PetscLogStageRegister("projectionStep");







