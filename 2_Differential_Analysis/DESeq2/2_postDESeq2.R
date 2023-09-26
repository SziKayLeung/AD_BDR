#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=50:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --error=test.e
#SBATCH --output=test.o

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
LOGEN_ROOT="/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen"
source ${LOGEN_ROOT}/miscellaneous
export PATH=$PATH:${LOGEN_ROOT}/miscellaneous

##-------------------------------------------------------------------------
inputBed=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/D_ONT/5_cupcake/8_characterise/bed12Files/ontBDR_concat_counts_coloured.bed12
sigiso=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/01_figures_tables/OntMerged_DESeq2_sigiso.txt
subset_fasta_gtf.py ${inputBed} --bed -i ${sigiso} -o OntMergedDESeq2Sig