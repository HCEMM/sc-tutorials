# Learn to use Slurm

So you can use our supercomputer the best way.

<br>

## Index

* [1. Who is this tutorial intended for](#1-who-is-this-tutorial-intended-for)

* [2. Starting point](#2-starting-point)

* [3. What is Slurm](#3-what-is-slurm)

* [4. Slurm commands](#4-slurm-commands)

* [5. Pratical examples of submitting jobs](#5-pratical-examples-of-submitting-jobs)

<br>

## 1. Who is this tutorial intended for

Everyone that would like to used the Scientific Computing cluster. Additionally, if you need to use Slurm on some other cluster that is managed by Slurm, the instructions in this tutorial will probably also be of use.

If you have any suggestions or questions, please [raise an issue](https://github.com/HCEMM/scientific_computing_docs/issues/new/choose).

<br>

## 2. Starting point

Say you have a Python script, handling big datasets, and your system runs out of memory while handling it (if you want to load a 20 Gb file, you need at least 20 Gb of memory).

Say you want to analyse NGS datasets, but they take too much of your storage. And while doing read alignment, you also run out of memory, because the databases/indexes are huge.

Say you want to do a molecular dynamics simulation, but it takes many days for simulating something simple, and you want to upscale it dramatically. Also, you have no GPUs on your system, and you would like to use GPUs to speed up things as much as possible, and maybe also do a little of deep-learning as a light hobby.

The Scientific Computing cluster is an appropriate response for all these challenges. If you want to learn more about the system itself, we suggest you take a look at [what]() was presented on the 15th of January 2025. You can also read the relevant [documentation file]().

<br>

## 3. What is Slurm

[Slurm](https://slurm.schedmd.com/quickstart.html) is the software responsible for distributing the computational resources available to the requests submitted. 

Slurm requires extensive information on the architecture of the system/cluster, which must be provided by the administrator of the system through the [configuration file](https://slurm.schedmd.com/slurm.conf.html).

<br>

## 4. Slurm commands

Slurm commands can be roughly grouped into functions for:
* starting/altering computational tasks (i.e., jobs, tasks, triggers, files transfer): `salloc`, `sattach`, `sbatch`, `sbcast`, `scancel`, `srun`, `strigger`
* supervision of resources usage and current state: `sacct`, `scontrol`, `sinfo`, `sprio`, `sshare`, `sstat`, `sview`

<br>

### 4.1. sbatch

To run an R script through Slurm, consider the following files:

`simple_script.R`
```R
get_hostname <- function() { 
  return (Sys.info()["nodename"])
}

print(paste("I am in server", get_hostname()))
```

`simple_script.sh`
```sh
#!/bin/bash
#SBATCH --job-name=multi-step-gpu-job      # Job name
#SBATCH --output=multi-step-gpu-job.out     # Output file for STDOUT
#SBATCH --error=multi-step-gpu-job.err      # Output file for STDERR
#SBATCH --ntasks=4                          # Number of tasks
#SBATCH --mem=16G                           # Memory per node (16 GB)
#SBATCH --time=02:00:00                     # Time limit (hh:mm:ss)
#SBATCH --partition=gpu                     # GPU partition to submit the job
#SBATCH --mail-type=ALL                     # Notification on job start, end, fail
#SBATCH --mail-user=your.email@example.com  # Where to send notifications

module load python/3.8                       # Load Python 3.8 for the environment

srun Rscript simple_script.R
``` 
After placing both files in either your HOME folder or your group's SCRATCH folder, you can run the R script through Slurm with the command `sbatch simple_script.sh`.

<br>

### 4.2. salloc

By running `salloc --account=jsequeira --partition=cpu --time=02:00:00`, you will submit a job request, requiring the default quantity of cores and memory, for two hours. When these resources become available, you will enter an interactive shell inside one of the compute nodes, and will be able to run commands directly on the compute node until you kill the job, or the time expires.

```console
jsequeira@master_node$ alloc --account=jsequeira --partition=cpu --time=02:00:00

salloc: Pending job allocation 92348
salloc: job 92348 queued and waiting for resources
salloc: job 92348 has been allocated resources
salloc: Granted job allocation 92348
salloc: Waiting for resource configuration
salloc: Nodes cpu5 are ready for job

jsequeira@cpu5$ 
```

See how the system changed from `master_node` to `cpu5`? This is a good way for users to test commands/scripts in the HPC environment, but shouldn't be used as the normal way of using the cluster. Ideally, you will prepare your scripts and submit them with the `sbatch` command, as shown before.

<br>

### 4.3. squeue

`squeue -u jsequeira` lists jobs waiting on the queue for user `jsequeira`.

`squeue -start` reports on the expected start time for pending jobs.

`squeue -j [jobID]` outputs the nodes executing a running job.

<br>

### 4.4. Other commands

`scancel [jobID]` cancels the job.

`sinfo` shows state of nodes and partitions.

`scontrol show [jobID]` shows detailed information about a job.

`sacct` displays accounting data for all jobs.

`sreport` shows accounting reports.

<br>

## 5. Pratical examples of Slurm

### 5.1. Tools installation

#### 5.1.1. Modules

The implementation of [modules](https://curc.readthedocs.io/en/latest/compute/modules.html) in Slurm provides a way of loading bundles of packages for specific tasks. **Module commands should not be run on the login node**, but rather in the batch scripts or interactive sessions - because these run in the compute nodes, where the modules are available. Here follows a list of some useful module commands:

* `module avail` lists all available modules.
* `module spider <module>` searches for a particular software.
* `module load <module>` loads a module into the system (e.g., `module load anaconda/2023.09`).
* `module unload <module>` removes the module.
* `module avail` lists all available modules.
* `module save <name>` saves the state of all loaded modules as collection `<name>`.
* `module restore <name>` restores the state of the modules to collection `<name>`.
* `module help` shows a man page for module commands.

#### 5.1.2. Spack

Has been designed for [reproducible builds](https://spack.readthedocs.io/en/latest/features.html) of software.

Can build [modules](https://spack-tutorial.readthedocs.io/en/latest/tutorial_modules.html).

#### 5.1.3. Conda environments

[Conda](https://docs.conda.io/en/latest/) is a great package manager that can provision packages for Python, C, and many other languages. By default, Conda provides a `base` environment, but additional environments can be created to better compartmentalize tools, and avoid conflicts with the base environment and the system itself.

```sh
#!/bin/bash
#SBATCH --job-name=conda-env-job            # Job name
#SBATCH --output=conda-env-job.out          # Standard output
#SBATCH --error=conda-env-job.err           # Standard error
#SBATCH --ntasks=1                          # Number of tasks
#SBATCH --cpus-per-task=2                   # CPUs per task
#SBATCH --mem=4G                            # Memory per task
#SBATCH --time=01:00:00                     # Time limit
#SBATCH --partition=cpu                     # Partition

# Load the Conda module
module load anaconda/2023.09

# Create a new Conda environment
conda create -n myenv python=3.8

# Activate the Conda environment
conda activate myenv

# Install some Python packages
conda install numpy pandas matplotlib

# Run the Python script
python my_script.py

# Deactivate the Conda environment
conda deactivate
```

#### 5.1.4. Python virtual environments

[Python virtual environments](https://docs.python.org/3/library/venv.html) can be used to create isolated environments for Python projects. They can be created and activated in Slurm batch scripts, as shown below.

```sh
#!/bin/bash
#SBATCH --job-name=venv-job                  # Job name
#SBATCH --output=venv-job.out                # Standard output
#SBATCH --error=venv-job.err                 # Standard error
#SBATCH --ntasks=1                           # Number of tasks
#SBATCH --cpus-per-task=2                    # CPUs per task
#SBATCH --mem=4G                             # Memory per task
#SBATCH --time=01:00:00                      # Time limit
#SBATCH --partition=cpu                      # Partition

# Create a new Python virtual environment
python -m venv myenv

# Activate the Python virtual environment
source myenv/bin/activate

# Install some Python packages
pip install numpy pandas matplotlib

# Run the Python script
python my_script.py

# Deactivate the Python virtual environment
deactivate
```

### 5.2. Running exclusive jobs

If running patient data, it might be important, in a legal sense, to ensure that data analysis is not exposed to other users. This can be achieved by running jobs in **exclusive** mode, where the job will be the only one running on the node.

For that, you simply need to add the `--exclusive` flag when running the `sbatch` command: 
```
sbatch --exclusive my_script.sh
```

In alternative, this flag can also be added to the batch script itself:

```sh
#SBATCH --exclusive
```

### 5.3. Multi-node batch scripts

This repository provides simple examples of Python scripts to run jobs whose **steps** can be spread out between **multiple nodes** of the same partition.  

#### 5.3.1. Job Array submission

[Job Array](https://slurm.schedmd.com/job_array.html) submission involves splitting input data into multiple parts, which are submitted for different runs of the same type of analysis. This approach might be very helpful when analysing either a very big dataset which can be analysed in separate parts, or analysing multiple datasets with the same approach.

Below is the example of a BATCH script for submitting a Job Array.

```shell
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
```

Slurm is able to split this job into several tasks (one for each run of `my_script.py`), which can then be run in parallel in multiple nodes.

#### 5.3.2. `&` signifies tasks as independent and parallelizable

When different datasets are to be submitted for different analyses/scripts, it can still be advantageous to submit them together in a single job. Parallelization, i.e., running the multiple tasks at the same time, will still be assured if connecting the independent tasks with `&`.

```sh
#!/bin/bash
#SBATCH --job-name=multi-task-parallel      # Job name
#SBATCH --output=multi-task-parallel.out    # Standard output
#SBATCH --error=multi-task-parallel.err     # Standard error
#SBATCH --ntasks=4                          # Number of parallel tasks
#SBATCH --cpus-per-task=2                   # CPUs per task
#SBATCH --mem=4G                            # Memory per task
#SBATCH --time=01:00:00                     # Time limit
#SBATCH --partition=cpu                     # Partition

# Load necessary modules
module load python/3.8

# Run each Python task in parallel using srun
srun --ntasks=1 --cpus-per-task=2 python my_script1.py --input input1.txt &
srun --ntasks=1 --cpus-per-task=2 python my_script2.py --input input2.txt &
srun --ntasks=1 --cpus-per-task=2 python my_script3.py --input input3.txt &
srun --ntasks=1 --cpus-per-task=2 python my_script4.py --input input4.txt &

# Wait for all tasks to complete
wait
```
Note the use of `wait`, which ensures the script waits for all the parallel tasks to finish before exiting.

#### 5.3.3. A multi-node job where tasks communicate

It is not trivial to run a job on Slurm capable of distributing multiple tasks into multiple servers, while keeping the tasks in communication. This can only be achieved with [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface).

What follows is a very simple example. Consider the scripts:

`mpi_parallel_script.py`

```py
from mpi4py import MPI

# Initialize the MPI communicator
comm = MPI.COMM_WORLD
rank = comm.Get_rank()   # Get the rank (task ID)
size = comm.Get_size()   # Get the total number of tasks

# Each process prints its rank and total size
print(f"Hello from task {rank} out of {size} tasks")

# Perform a simple computation: each task computes its square
data = rank ** 2
print(f"Task {rank} computed value: {data}")

# Gather the results from all tasks (root is rank 0)
results = comm.gather(data, root=0)

if rank == 0:
    print(f"Results gathered at root: {results}")
```

`mpi_parallel_script.sh`

```sh
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
```
By submitting this job with `sbatch mpi_parallel_script.sh`, you will be running a multi-node job with 4 tasks, where each task will print its rank and total size, compute its square, and gather the results at the root task (rank 0). 

<br>




