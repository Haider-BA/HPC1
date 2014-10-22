#include "headers.h"

void printvec(Vec X, const char file[])
{
	PetscViewer me;
	PetscViewerCreate(PETSC_COMM_WORLD, &me);
	PetscViewerFileSetName(me, file);
	PetscViewerSetFormat(me, PETSC_VIEWER_ASCII_SYMMODU);
	VecView(X, me);
	PetscViewerDestroy(me);
}

