#include "../NavierStokesSolver.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

PetscErrorCode NavierStokesSolver::writeData()
{
	PetscErrorCode ierr;
	PetscInt	rank;
	Vec		qxGlobal, qyGlobal, phi;
//	PetscReal	**uGlobal;
	PetscViewer	viewer;

	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if(rank==0)
	{
		PetscInt its1, its2;
		std::string filename = "cases/iterationCount.txt";
		std::ofstream	iterationFile;
		if(timeStep==1)
		{
			iterationFile.open(filename.c_str());
		}
		else
		{
			iterationFile.open(filename.c_str(), std::ios::out | std::ios::app);
		}
		ierr = KSPGetIterationNumber(ksp1,&its1);
		ierr = KSPGetIterationNumber(ksp2,&its2);
		iterationFile << timeStep << '\t'<< its1 << "\t" << its2 << std::endl;
		iterationFile.close();
	
		std::ofstream f("cases/simulationInfo.txt");
		f << "-nx\t" << fluid.nx << '\n';
		f << "-ny\t" << fluid.ny << '\n';
		f << "-startStep\t"<< 0 << '\n';
		f << "-nt\t" << fluid.nt << '\n';
		f << "-nsave\t" << fluid.nsave << '\n';
		f << "-dt\t" << fluid.dt << '\n';
		f.close();
	}

     if(timeStep%fluid.nsave ==0)
     {
	std::string fileName, savePointDir;
	std::stringstream ss;
	ss << "cases/" << std::setfill('0') << std::setw(7) << timeStep;
	savePointDir = ss.str();
	mkdir(savePointDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);	

	ierr = DMCompositeGetAccess(qPack, q, &qxGlobal, &qyGlobal);

	ss.str("");
	ss.clear(); // ?
	ss << savePointDir << "/qx.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,  fileName.c_str(), FILE_MODE_WRITE, &viewer);
	ierr = VecView(qxGlobal, viewer);
//	ierr = DAVecGetArray(uda,qxGlobal, &uGlobal);
//	ierr = PetscPrintf(PETSC_COMM_WORLD,"qx vector"); 
//	ierr = VecView(qxGlobal, PETSC_VIEWER_STDOUT_WORLD);	
	ierr = PetscViewerDestroy(&viewer);

	ss.str("");
	ss.clear();
	ss << savePointDir << "/qy.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_WRITE, &viewer);
	ierr = VecView(qyGlobal, viewer);
//	ierr = PetscPrintf(PETSC_COMM_WORLD,"qy vector");
//	ierr = VecView(qyGlobal,PETSC_VIEWER_STDOUT_WORLD);
	ierr = PetscViewerDestroy(&viewer);
	ierr = DMCompositeRestoreAccess(qPack, q, &qxGlobal, &qyGlobal);
	
	
	ss.str("");
	ss.clear();
	ss << savePointDir << "/phi.dat";
	fileName = ss.str();
	ierr = DMCompositeGetAccess(lambdaPack, lambda, &phi);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_WRITE, &viewer);
	ierr = VecView(phi, viewer);
	ierr = PetscViewerDestroy(&viewer);
	
	ierr = DMCompositeRestoreAccess(lambdaPack, lambda, &phi);

	ierr = PetscPrintf(PETSC_COMM_WORLD, "data written done\n");
     }

     return 0;
}
