#include <iostream>
#include <cmath>

using namespace std;

int main()
{
	const int N= 256*256;	
	int nx= 256;
	int ny = 256;
	double h = 1./(nx-1);
	int MaxStep = int(nx*nx*10/M_PI/M_PI) + 100;
	double convergence_epsilon = 10e-5;

	double new_u[N] = {};
	double initial_u[N] = {};
	double ext_u[N] = {};

	for(int I=0; I<N; I++)
	{
		double v=0.0;
		int i=I/nx;    // y=i*h
		int j=I - i*nx; // x=j*h

		if(i==0) v= v + sin(M_PI*j*h);
		if(i==nx-1) v=v+sin(M_PI*j*h) * exp(-M_PI);
		if( v!=0) { initial_u[I] = v; new_u[I] = v; };
		v= sin(M_PI*j*h) * exp(-M_PI*i*h);
		ext_u[I] = v;
	}

	double error[N] = {}; 
	for(int istep=0; istep<MaxStep; istep++)
	{
		double ErrorSum = 0.0;
		for(int i=1; i<nx-1; i++)
		{
			for(int j=1; j<ny-1; j++)
			{
				int I = i*nx + j; 
				new_u[I] = (initial_u[I-1] + initial_u[I+1] + initial_u[I+nx] + initial_u[I-nx])/4;
			

		error[I] = new_u[I] - initial_u[I];
		ErrorSum += error[I] * error[I];
			}
		}

	ErrorSum = sqrt(ErrorSum);
	cout << "step " << istep << " " << "error " << ErrorSum << endl;
		
	if(ErrorSum < convergence_epsilon)
	{
		cout << istep << endl;
		std::cout << "Jacobi Iteration Convergence " << std::endl;
		break;
	}
	else
	{
		if(istep==MaxStep)
		{
			std::cout << "Max Iteration Steps have reached, please increase Steps " << std::endl;
		}
		for(int I=0; I<N; I++)
		{
			initial_u[I] = new_u[I];
		}
	}
	}

}

