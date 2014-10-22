#include "headers.h"

Params parameters(){
	
	Params P;
	P.Nx = 100;
	P.Ny = 102;

	P.hx = 1.00/(P.Nx);
	P.hy = 1.00/(P.Ny - 1);

	return P;
}

//by default construct,  Params P = parameters();


