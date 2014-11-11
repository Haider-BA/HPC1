#include "JacobianSolver.h"
#include <cmath>

#ifdef	_ViewVect
	#include <iostream>
	using namespace std;
#endif


void JacobianSolver::VectAssembly()
{
/*
	rhs = new double[N];
	extsol = new double[N];
*/
	for(int I=0; I<N;I++)
	{
		double v=0.0;
		int i=I/nx; // y = i*h
		int j=I - i*nx; // x = j*h

		if(i==0) v = v + sin(M_PI*j*h);
		if(i==nx-1) v = v + sin(M_PI*j*h) * exp(-M_PI);
		if(j==0) v = 0.0;
		if(j==ny-1) v = 0.0;
		
		if(v != 0.0) { rhs[I] = v; };
		v = sin(M_PI*j*h) * exp(-M_PI*i*h);
		extsol[I] = v;
	}

#ifdef _ViewVect
	for(int i=0; i<N; i++)
	{
		cout << rhs[i] << " " ;
	}
	cout << endl;
	
	for(int i=0; i<N; i++)
	{
		cout << extsol[i] << " ";
	
	}
	cout << endl;
#endif

}
