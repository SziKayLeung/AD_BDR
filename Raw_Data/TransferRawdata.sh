#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=60:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --reservation=research_project-mrc148213_5

# 08/05/2021: Transfer rawdata from lims to ISCA 

cd /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Rawdata
wget -m ftp://Project_10390:qM9rgx6A9yw3@ftp1.sequencing.exeter.ac.uk/
mv ftp1.sequencing.exeter.ac.uk/* .
rmdir ftp1.sequencing.exeter.ac.uk
