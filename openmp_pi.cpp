#include<iostream>  
#include<iomanip>
#include<cmath> 
#include<time.h>
#include<omp.h>  

using namespace std;

	int main()
{
	double	minus_one = -1;

	for(long int N=10e4; N < 1e9; N *= 10)
	{	
		int i;
		double sum = 0;		
		clock_t t1=clock();
		#pragma omp parallel for private(i) reduction(+:sum) 
		for(i=0; i<N; i++){
		sum +=  pow(minus_one, i)/(2*i + 1);
		}	
		clock_t t2=clock();
		double cost = (double) (t2 -t1)/1e6;
		cout << "when N = " << N << ", Pi = "<<  setprecision(10) << 4*sum << " Running Time = " << cost <<endl;
	}

	return 0;
}


