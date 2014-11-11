#ifndef JacobianSolver_h
#define JacobianSolver_h

class JacobianSolver
{
	private:
		const int nx, ny;
		const int N;
		int MaxSteps;
		const double h; 
		double** Dinverse;
		double** LU;
		double* rhs;
		double* sol;
		double* lastsol;
	        double* extsol;
		double convergence_epsilon;
	
	public:
	JacobianSolver(int nx, int ny, double h)
		
		: nx(nx),
		  ny(ny),
		  N(nx*ny),
		  h(h)
		  {
			  Dinverse = new double* [N];
			  for(int i=0; i<N; i++)
			  {
				  Dinverse[i] = new double[N];
			  }

			  LU = new double* [N];
			  for(int i=0; i<N; i++)
			  {
				  LU[i] = new double[N];
			  }

			  rhs = new double[N];
			  extsol = new double[N];

//			  sol = new double[N];
//			  lastsol = new double[N];

		  }

	~JacobianSolver()
	{
		for(int i=0; i<N; i++)
		{
			delete[] Dinverse[i];
			delete[] LU[i];
		}
		delete[] Dinverse;
		delete[] LU;

		delete[] rhs;
		delete[] extsol;
//		delete[] sol;
//		delete[] lastsol;
	}
		
	void MatAssembly();
	void VectAssembly();
	void Jacobian();
};

#endif
