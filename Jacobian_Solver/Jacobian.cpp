#include "JacobianSolver.h"
#include <iostream>
#include <cmath>

using namespace std;

void JacobianSolver::Jacobian()
{
	MaxSteps=10;
	convergence_epsilon=10e-5;

      	for(int istep=0; istep<MaxSteps; istep++)
  {
	sol = new double[N];
	lastsol = new double[N];
	double* Error = new double[N];
	if(istep==MaxSteps) std::cout << "Reach Max Iterations " << std::endl; 
  	double L2norm(0.0f);
	for (int i=0; i< nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			sol[i] = Dinverse[i][j]*(LU[i][j]*lastsol[i] + rhs[i]);
		}
		Error[i]  = sol[i] - extsol[i];
//		Error[i] = sol[i] - lastsol[i];	
		lastsol[i] = sol[i]; 
		L2norm += Error[i] * Error[i];
	}
	
	cout << " ----------------- " << endl;

	cout << "ext vector " << endl;
	for(int i=0; i<N; i++)
	{
		cout << extsol[i] << " ";
	}
	cout << endl;	
		

	cout << "sol vector " << endl;
	for(int i=0; i<N; i++)
	{
		cout << sol[i] << " ";
	}
	cout << endl;	
		
	cout<< "Error " << endl;
	for(int i=0; i<N; i++)
	{
		cout << Error[i] << " ";
	}
	cout << endl;
	

		L2norm = sqrt(L2norm);	
		std::cout<< "step " << istep << " L2norm " << L2norm << std::endl;


	if( L2norm < convergence_epsilon)
	{
		std::cout << "Jacobian Converged " << std::endl;
	}

	delete[] Error;
	delete[] lastsol;
	delete[] sol;
   }
}


		



