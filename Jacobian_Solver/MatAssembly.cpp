#include "JacobianSolver.h"

#ifdef	_ViewMat
	#include <iostream>
	using namespace std;
#endif

void JacobianSolver::MatAssembly()
{
	double v;
	int J; 
	for( int I=0; I<N; I++)
	{
		v=-1.0;
		int i = I/nx;
		int j = I%nx;

		if(i>0) {J=I-nx; LU[I][J] = v;}
		if(i<ny-1) {J=I+nx; LU[I][J] = v;}
		if(j>0) { J=I-1; LU[I][J] = v;}
		if(j<nx-1) {J=I+1; LU[I][J]=v;}

		v=1./4.0;
		Dinverse[I][I]=v;
	}

#ifdef _ViewMat
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<N; j++)
		{
			cout<< LU[i][j]+ 1.0/Dinverse[i][j] << " ";
		}
		cout << endl;
	}
#endif

}



