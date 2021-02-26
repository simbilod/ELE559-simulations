#!/bin/bash
#SBATCH --job-name=ECE559-bent   # create a short name for your job
#SBATCH --output=slurm-%A.%a.out # STDOUT file
#SBATCH --error=slurm-%A.%a.err  # STDERR file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=4        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=00:10:00          # total run time limit (HH:MM:SS)
#SBATCH --array=0-2              # job array with index values 0, 1, 2
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=<your-net-id>@princeton.edu

# Read more about slurm at https://researchcomputing.princeton.edu/slurm

source /home/ELE559/ELE559.bashrc

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

mpirun -np 4 python bent_waveguide.py