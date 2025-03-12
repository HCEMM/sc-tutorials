#!/bin/bash
#SBATCH --job-name=multi-step-gpu-job      # Job name
#SBATCH --output=multi-step-gpu-job.out     # Output file for STDOUT
#SBATCH --error=multi-step-gpu-job.err      # Output file for STDERR
#SBATCH --ntasks=4                          # Number of tasks
#SBATCH --gpus=2                            # Request 2 GPUs per node
#SBATCH --nodes=1                           # Number of nodes (1 in this case)
#SBATCH --cpus-per-task=4                   # Number of CPUs per task
#SBATCH --mem=16G                           # Memory per node (16 GB)
#SBATCH --time=02:00:00                     # Time limit (hh:mm:ss)
#SBATCH --partition=gpu                     # GPU partition to submit the job
#SBATCH --mail-type=ALL                     # Notification on job start, end, fail
#SBATCH --mail-user=your.email@example.com  # Where to send notifications

# Load necessary modules
module load cuda/11.2                        # Load CUDA module (version 11.2)
module load python/3.8                       # Load Python 3.8 for the environment

# Define different job steps with specific tasks

# Step 1: Data preprocessing (CPU-only task)
srun --ntasks=1 --cpus-per-task=2 --mem=4G \
     --job-name="data-preprocessing" \
     python preprocess_data.py --input data/raw --output data/processed

# Step 2: Model training (GPU-accelerated task)
srun --ntasks=1 --gpus=2 --cpus-per-task=4 --mem=8G \
     --job-name="model-training" \
     python train_model.py --data data/processed --epochs 50 --batch-size 32

# Step 3: Model evaluation (GPU-accelerated task)
srun --ntasks=1 --gpus=1 --cpus-per-task=4 --mem=4G \
     --job-name="model-evaluation" \
     python evaluate_model.py --data data/processed --model output/model.pt

# Step 4: Postprocessing results (CPU-only task)
srun --ntasks=1 --cpus-per-task=2 --mem=2G \
     --job-name="postprocessing-results" \
     python postprocess_results.py --input output/results --output final_report.csv

echo "Job complete!"
