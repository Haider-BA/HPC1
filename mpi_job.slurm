#!/bin/sh
#SBATCH --partition=general-compute
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --constraint=CPU-E7-4830
#SBATCH --mem=48000
#SBATCH --job-name="mpi_pi"
#SBATCH --output=mpi_pi.out
#SBATCH --mail-user=zhengjia@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --requeue

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST="$SLURM_JOB_NODELIST
echo "SLURM_NNODES="$SLURM_NNODES
module load intel
module load intel-mpi
ulimit -s unlimited

for ((np=1; np<=32; np *= 2)) do
	echo "runing mpi_pi on $np cpus"
	mpiicpc -o mpi_pi mpi_pi.cpp
	mpirun -n $np ./mpi_pi
done

