#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=144:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

##-------------------------------------------------------------------------

SC=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR/1_ONT_Pipeline

# command line - create folders
#bash ${SC}/0_create_folder_batch2.sh

# merge fastq5 pass in Batch2 
#jid1=$(sbatch ${SC}/1_merge_batch2.sh)
jid1=$(sbatch ${SC}/check.sh)

# batch and align samples separately
jid2=$(sbatch --dependency=afterok:$jid1 --job-name=demux ${SC}/2_demux_batch2.sh)

# merge all samples, collapse and run sqanti
jid3=$(sbatch --dependency=afterok:$jid2 --job-name=tclean ${SC}/3_cutadapt_minimap2_tclean_batch2.sh)
