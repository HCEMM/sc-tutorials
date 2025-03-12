#!/bin/bash
#SBATCH --job-name=parallel-tasks           # Job name
#SBATCH --output=output_%A_%a.out           # Standard output (A is the job ID, a is the array index)
#SBATCH --error=error_%A_%a.err             # Standard error
#SBATCH --array=0-3                         # Job array with 4 tasks (index 0, 1, 2, 3)
#SBATCH --ntasks=1                          # One task per array element
#SBATCH --cpus-per-task=2                   # Number of CPUs per task
#SBATCH --mem=4G                            # Memory per task
#SBATCH --time=01:00:00                     # Time limit
#SBATCH --partition=cpu                     # Partition

# Load necessary modules
module load python/3.8

# Define different inputs or parameters for each array element
input_files=("input1.txt" "input2.txt" "input3.txt" "input4.txt")

# Run the Python script with the input corresponding to the array task ID
python my_script.py --input ${input_files[$SLURM_ARRAY_TASK_ID]}