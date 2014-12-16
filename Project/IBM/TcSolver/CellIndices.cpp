#include "../TcSolver.h"

PetscErrorCode TcSolver::calculateCellIndices()
{
	std::vector<PetscReal> &xMesh = NavierStokesSolver::fluid.x,
			       &yMesh = NavierStokesSolver::fluid.y;

	I.reserve(x.size()); // x is for body
	J.reserve(x.size());

	PetscInt i=0, j=0;

	//find cell for zeroth point
	while(xMesh[i+1] < x[0]) i++;
	while(yMesh[j+1] < y[0]) j++;

	I.push_back(i);
	J.push_back(j);

	for(size_t l=1; l<x.size(); l++)
	{
		//if the next boundary point is to the left of the current boundary point
		if(x[l] < x[l-1])
		{
			while(xMesh[i] > x[l])
				i--;  // while ture, means the next right grid-node is also inside the body
		}
		//if the next boundary point is to the right of the current boundary point
		else
		{
			while(xMesh[i+1] < x[l])
				i++;
		}	

// x/yMesh belongs to fluid; x/y belongs to body
		//if the next boundary point is below the current boundary point
		if(y[l] < y[l-1])
		{
			while(yMesh[j] > y[l])
				j--;
		}
		//if the next boundary point is above the current boundary point
		else
		{
			while(yMesh[j+1] < y[l])
				j++;
		}

		I.push_back(i);
		J.push_back(j);

	}
	return 0;
}
