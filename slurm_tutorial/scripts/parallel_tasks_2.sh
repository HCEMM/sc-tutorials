#!/bin/bash
#SBATCH --job-name=multi-task-parallel     # Job name
#SBATCH --output=multi-task-parallel.out   # Standard output
#SBATCH --error=multi-task-parallel.err    # Standard error
#SBATCH --ntasks=4                         # Number of parallel tasks
#SBATCH --cpus-per-task=2                  # CPUs per task
#SBATCH --mem=4G                           # Memory per task
#SBATCH --time=01:00:00                    # Time limit
#SBATCH --partition=standard               # Partition

# Load necessary modules
module load python/3.8

# Run each Python task in parallel using srun
srun --ntasks=1 --cpus-per-task=2 python my_script.py --input input1.txt &
srun --ntasks=1 --cpus-per-task=2 python my_script.py --input input2.txt &
srun --ntasks=1 --cpus-per-task=2 python my_script.py --input input3.txt &
srun --ntasks=1 --cpus-per-task=2 python my_script.py --input input4.txt &

# Wait for all tasks to complete
wait
