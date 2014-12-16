#include "../TcSolver.h"

PetscBool TcSolver::isInfluenced(PetscReal xGrid, PetscReal yGrid, PetscReal xBody, PetscReal yBody, PetscReal radius, PetscReal *disp)
{
	PetscReal width[2];
	PetscReal nx = NavierStokesSolver::fluid.nx,
		  ny = NavierStokesSolver::fluid.ny;

	std::vector<PetscReal> &xMesh = NavierStokesSolver::fluid.x,
			       &yMesh = NavierStokesSolver::fluid.y;

	width[0] = xMesh[nx]-xMesh[0];
	width[1] = yMesh[ny]-yMesh[0];

	disp[0] = fabs(xGrid - xBody);
	disp[1] = fabs(yGrid - yBody);

	return (disp[0] < radius && disp[1] < radius) ? PETSC_TRUE : PETSC_FALSE;

}
