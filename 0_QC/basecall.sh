#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=S.K.Leung@exeter.ac.uk # email address
#SBATCH --output=basecall.o
#SBATCH --error=basecall.e

module load Miniconda2 
source activate sqanti2_py3

cd /lustre/projects/Research_Project-MRC148213/sl693/AD_BDR/1_raw/Batch2/
pod5 convert fast5 ./fast5_pass/*.fast5 --output pod5_pass --one-to-one ./fast5_pass/