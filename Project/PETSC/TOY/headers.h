#include<stdlib.h>
#include "petscsnes.h"
#include "petscda.h"
#include "petsctime.h"

typedef struct{
	PetscScalar u, v;
} Field;

// parameters in the problem, value be set in .c
typedef struct{
	PetscInt Nx, Ny;
	PetscScalar hx, hy;
} Params;

//data passed to functions when evaluation equations and Jacobian
typedef struct{
	Params P;
	DA da;
} Data;

//functions in parameter.c
Params parameters();

//functions in utilities.c
void printvec(Vec X, const char file[]);


