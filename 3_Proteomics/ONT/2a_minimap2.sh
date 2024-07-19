#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrchq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=2a_minimap2.o
#SBATCH --error=2a_minimap2.e

# 10/10/2023: Minimap alignment of best orf to genome fasta for isoform visulisation

#-----------------------------------------------------------------------#

# load and source packages
module load Miniconda2

source activate lrp
source /lustre/projects/Research_Project-MRC148213/lsl693/scripts/General/4_Proteogenomics/1_proteogenomics_functions.sh
source /lustre/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR/3_Proteomics/2a_bdr_proteomics.config

# run minimap2
source activate nanopore  
cd $WKD_ROOT/5_calledOrfs
minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 $GENOME_FASTA $NAME"_bestORF.fasta" > $NAME.sam 2> $NAME.map.log
samtools sort -O SAM $NAME.sam > $NAME.sorted.sam