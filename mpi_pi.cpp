#include<iostream>
#include<iomanip>
#include<cmath>
#include<time.h>
#include<mpi.h>		

using namespace std;

int main()
{
	MPI_Init(NULL,NULL);
	
	for(long int N = 10e4; N<1e9; N *=10){
	int num_proc = 32; 	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	double global_sum = 0.0;
	double local_sum = 0.0;


	int num_elements_per_proc = int(N/num_proc);
	
	int local_min_index = world_rank * num_elements_per_proc;
	int local_max_index = (world_rank+1) * num_elements_per_proc;
	
	int minus_one = -1;

	clock_t t1 = clock();

	for (int i=local_min_index; i < local_max_index; i++)
	{
		local_sum += pow(minus_one, i)/(2*i + 1);
	}		
	
	MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	if(world_rank ==0){
	cout << "pi = " << setprecision(10) << global_sum*4 << endl;
	clock_t t2 = clock();
	double cost = (double) (t2 - t1)/1e6;
	cout << "when N= " << N << " and num of proc = " << num_proc <<  "running time of MPI " << cost << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}
	MPI_Finalize();

		return 0;
}
 	
