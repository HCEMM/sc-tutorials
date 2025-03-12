#!/bin/bash
#SBATCH --job-name=mpi-multi-task-job       # Job name
#SBATCH --output=mpi-output-%j.out          # Standard output (one file for all tasks)
#SBATCH --error=mpi-error-%j.err            # Standard error (one file for all tasks)
#SBATCH --ntasks=4                          # Request 4 tasks
#SBATCH --nodes=2                           # Use 2 nodes (can span nodes)
#SBATCH --cpus-per-task=2                   # 2 CPUs per task
#SBATCH --mem-per-cpu=2G                    # Memory per CPU (2 GB per CPU)
#SBATCH --time=01:00:00                     # Time limit (1 hour)
#SBATCH --partition=standard                # Partition to submit to

# Load necessary modules
module load python/3.8                      # Python module
module load mpi/openmpi-4.0.5               # MPI library (OpenMPI in this case)

# Run the Python script using MPI
mpirun -n 4 python mpi_parallel_script.py
