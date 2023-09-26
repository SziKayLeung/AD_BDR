#!/bin/bash

## make folder and move previous results from Batch1 to a designated folder

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh

# previous run 
cd ${WKD_ROOT}
mkdir -p 1_demultiplex/Batch1 1_demultiplex/Batch2
mv 1_demultiplex/*PAM95037* 1_demultiplex/Batch1

mkdir -p 1b_demultiplex_merged/Batch1 1b_demultiplex_merged/Batch2 1b_demultiplex_merged/Batch2/log
mv 1b_demultiplex_merged/*BC* 1b_demultiplex_merged/Batch1 
mv 1b_demultiplex_merged/log 1b_demultiplex_merged/Batch1 

mkdir -p 2_cutadapt_merge/Batch1 2_cutadapt_merge/Batch2
mv 2_cutadapt_merge/*BC* 2_cutadapt_merge/Batch1

mkdir -p 3_minimap/Batch1 3_minimap/Batch2
mv 3_minimap/*BC* 3_minimap/Batch1

mkdir -p 4_tclean/Batch1 4_tclean/Batch2
mv 4_tclean/*BC* 4_tclean/Batch1

mkdir 5_talon
mv 5_talon_label 5_talon/
mv 6_talon 5_talon/
mv 7_sqanti3/ 5_talon/

mkdir 5_cupcake 5_cupcake/5_align
mkdir -p 5_cupcake/5_align/Batch1 5_cupcake/5_align/Batch2 