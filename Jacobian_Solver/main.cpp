#include "JacobianSolver.h"


int main()
{
	double domainx(1.0), domainy(1.0);
//        double h0[] ={ 1./3., 1./6., 1./12.};
//	for(int i=0; i<sizeof(h0)/sizeof(h0[0]); i++)

	{
//		double h=h0[i];
		double h = 1./4.;
		int nx = (int) domainx/h;
		int ny = (int) domainy/h;
		JacobianSolver hw4(nx, ny, h);
		
		hw4.MatAssembly();
		hw4.VectAssembly();
		hw4.Jacobian();
	}

	return 0;
}


