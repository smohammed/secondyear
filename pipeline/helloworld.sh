#!/bin/sh
#
# Simple "Hello World" submit script for Slurm.
#
# Replace <ACCOUNT> with your account name before submitting.
#
#SBATCH --account=<ACCOUNT>      # The account name for the job.
#SBATCH --job-name=HelloWorld    # The job name.
#SBATCH -c 1                     # The number of cpu cores to use.
#SBATCH --time=1:00              # The time the job will take to run (here, 1 min)
#SBATCH --mem-per-cpu=1gb        # The memory the job will use per cpu core.
 
echo "Hello World"
sleep 10
date
 
# End of script