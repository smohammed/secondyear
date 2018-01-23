#!/bin/sh
#
# Simple "Hello World" submit script for Slurm.
#
#SBATCH --account=astro # The account name for the job.
#SBATCH --job-name=GALEX # The job name.
#SBATCH -c 4 # The number of cpu cores to use.
#SBATCH --time=60:00 # The time the job will take to run.
#SBATCH --mem-per-cpu=5gb # The memory the job will use per cpu core.
 
module load anaconda
 
#Command to execute Python program
python smoothscans_cluster.py
 
#End of script