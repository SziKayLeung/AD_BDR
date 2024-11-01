#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=3:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

## ---------------------------
## Purpose: Generate files (Iso-Seq annotation and expression) needed for tappAS 
## 
## 05/08/2022: Recreated files
## ---------------------------


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
export SC_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
source $SC_ROOT/1_IsoSeq_Pipeline/bdr_isoseq.config
source $SC_ROOT/2_Differential_Analysis/bdr_differential.config
source $SC_ROOT/2_Differential_Analysis/01_source_functions.sh

##-------------------------------------------------------------------------
#************************************* Prepare files for TappAS
# 4 files:
# 1) Iso-Seq Annotation scaffold file                       = XX.collapsed.gff3
# 2) Iso-Seq targeted isoforms                              = XX_TargetTrans.txt
# 3) Iso-Seq Expression File
# 4) Iso-Seq Phenotype File


# create directory 
mkdir -p $ISO_DIFF_DIR

# File 1
# Annotation file generated from IsoAnnotLite in SQANTI3
# Note, RNA-Seq reads were not used as junction filter, given lower depth and reduce false negative
cp $ISO_WKD_ROOT/$NAME".collapsed.gff3" $TAPPAS_INPUT_DIR

# File 2
# Tab file of arget isoforms associated to target gene to be used as inclusion criteria for tappAS
# This means that counts from non-target isofoms will still be used for normalisation
# but only the relevant target isoforms will be used for downstream purposes (thereby speeding tappAS)
subset_targets $ISO_DIFF_DIR

# File 3
# regenerate expression matrix for tappas from long reads FL read counts, 
# include non-target isoforms to ensure not skewing downstream normalisation
# counts_subset_4tappas <input_class> <output_class> <type_genes>
counts_subset_4tappas $ISO_WKD_ROOT/$SQNAME"_classification.txt" $ISO_DIFF_DIR/$NAME".expression.txt" notAD

# File 4
cp $ISO_TAPPAS_PHENO $ISO_DIFF_DIR

##-------------------------------------------------------------------------
#************************************* Transfer files 
mkdir -p $TAPPAS_INPUT_DIR/IsoSeq_Expression
cp $ISO_DIFF_DIR/* $TAPPAS_INPUT_DIR/IsoSeq_Expression
