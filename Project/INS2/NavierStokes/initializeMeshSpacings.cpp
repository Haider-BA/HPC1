// for staggered mesh, cell for x-momentom, y-momentom is different
//  so here make differnt vector to store the x/y-spacing for each momentum

#include "../NavierStokesSolver.h"

void NavierStokesSolver::initializeMeshSpacings()
{
	PetscInt numX, numY;

	//mesh spacing for U
	numX = fluid.nx - 1; 
	numY = fluid.ny;

	//dx
	dxU.resize(numX+1); // here including the spacing to XMinus/XPlus boundary

	for(PetscInt i=0; i<numX+1; i++)
	{
		dxU[i] = fluid.dx;
	}
	//dy
	dyU.resize(numY+1); //consider two ghost b.c. nodes outside of the domain
	for(PetscInt j=1; j<numY;j++)
	{
		//check if the point is at top/bottom edge
		dyU[j] = fluid.dy;
	}
	dyU[0] = 0.5*fluid.dy;
	dyU[numY] = 0.5*fluid.dy;

	//mesh spacing for V
	numX = fluid.nx;
	numY = fluid.ny - 1;
	//dx
	dxV.resize(numX+1);
	for(PetscInt i=1; i<numX; i++)
	{
		dxV[i] = fluid.dx;
	}
	dxV[0] = 0.5*fluid.dx;
	dxV[numX]=0.5*fluid.dx;

	//dy
	dyV.resize(numY+1);
	for(PetscInt j=0; j<numY+1; j++)
	{
		dyV[j] = fluid.dy;
	}
}


