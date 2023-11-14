#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrchq # submit to the parallel queue
#SBATCH --time=100:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=2_run_proteogenomics.o
#SBATCH --error=2_run_proteogenomics.e
#SBATCH --mem=500G 

# 02/08/2022: Run proteogenomics pipeline on all 96 BDR samples

#-----------------------------------------------------------------------#
## print start date and time
echo Job started on:
date -u

module load Miniconda2
source activate lrp
source /gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen/proteomics/proteogenomics.sh
source /lustre/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR/3_Proteomics/2a_bdr_proteomics.config

echo "#************************************* Collate and prepare long-read data"
#collate_longread_processed
#prepare_reference_tables
#summarise_longread_data

echo "#************************************* Call open reading frames and classify proteins"
#call_orf
#refine_calledorf
#classify_protein

echo "#************************************* Run hybrid annotations"
#run_hybrid_annotation

echo "#************************************* Run Metamorpheus"
#run_metamorpheus filtered
#run_metamorpheus refined
#run_metamorpheus hybrid
#run_metamorpheus gencode
#run_metamorpheus uniprot
#run_metamorpheus_rescue

echo "#************************************* Generate output tracks"
run_peptide_analysis
generate_cds_tracks
generate_peptide_tracks

echo "#************************************* Other output data"
compare_protein_groups
identify_novel_peptides

echo "#***************All done!****************#"
