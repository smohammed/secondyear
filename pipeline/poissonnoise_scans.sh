#!/bin/sh
#
# Adding poisson noise to images
#
#SBATCH --account=astro # The account name for the job.
#SBATCH --job-name=GALEX # The job name.
#SBATCH --exclusive
#SBATCH -N 1
#SBATCH --time=10:00:00 # The time the job will take to run.
 
module load anaconda/2-4.2.0

#Command to execute Python program
python poissonnoise_scans.py

#End of script
