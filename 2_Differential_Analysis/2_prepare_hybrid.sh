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
# 2) RNA-Seq Expression File
# 3) RNA-Seq Phenotype File


# create directory 
mkdir -p $HYBRID_DIFF_DIR
mkdir -p $TAPPAS_INPUT_DIR/RNASeq_Expression

# File 1
# Same annotation file specified in 1_prepare_isoseq_solo.sh
cp $ISO_WKD_ROOT/$NAME".collapsed.gff3" $TAPPAS_INPUT_DIR

# File 2
# Align RNA-Seq files (n = 96) to Iso-Seq annotation
bash $SC/2b_align_rnaseq2isoseq.sh
wait

cp $HYBRID_DIFF_DIR/$NAME"_RNASeq.expression.txt" $TAPPAS_INPUT_DIR/RNASeq_Expression

# File 3
cp $HYBRID_TAPPAS_PHENO $TAPPAS_INPUT_DIR/RNASeq_Expression

##-------------------------------------------------------------------------
#************************************* Account for missing samples (run 2b_validate_rnaseq_samples.Rmd)
Rscript -e "rmarkdown::render('$SC_ROOT/2_Differential_Analysis/2c_validate_rnaseq_samples.Rmd')"
# modified expression file
cp $HYBRID_DIFF_DIR/$NAME"_RNASeq_Final.expression.txt" $TAPPAS_INPUT_DIR/RNASeq_Expression
# modified phenotype file
cp $TAPPAS_PHENO_DIR/ADBDR_MinusIntermediate_RNASeqPhenotypeTAPPAS.txt $TAPPAS_INPUT_DIR/RNASeq_Expression