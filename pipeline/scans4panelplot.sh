#!/bin/sh
#
# Simple "Hello World" submit script for Slurm.
#
#SBATCH --account=astro # The account name for the job.
#SBATCH --job-name=GALEX # The job name.
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --time=3:00:00 # The time the job will take to run.
 
module load anaconda
 
#Command to execute Python program
python scans4panelplot.py

#End of script