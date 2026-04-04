#!/bin/bash
#SBATCH --job-name=AFT                          # Job name
#SBATCH --array=1-1000                               # At most 50 jobs at a time YMMV, make sure 
#SBATCH --cpus-per-task=1
#SBATCH --account=statistics
#SBATCH --mail-type=ALL                            # Mail events (NONE, BEGIN, END, FAIL, AL$
#SBATCH --mail-user=woojung.bae@ufl.edu    # Where to send mail
#SBATCH --mem=10G                                   # Memory limit
#SBATCH --partition=hpg-default
#SBATCH --time=2-00:00:00                          # Time limit hrs:min:sec
#SBATCH --output=outfiles/array_%A-%a_AFT.log  # the SLURM outputs will be created i$

date; hostname; pwd

module load gcc R

#Run R script
Rscript package_sim2_hpc_TIME.R
