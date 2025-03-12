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
