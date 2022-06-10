#!/bin/bash
#SBATCH -N 1
#SBATCH -c 2
#SBATCH -p genacc_q
#SBATCH -t 60:00:00
#SBATCH --mem-per-cpu 3900M
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sj14m@my.fsu.edu


module load R/3.5.2



R CMD BATCH --no-save --no-restore  2nd_revised_datachecking_CTG_new.R
