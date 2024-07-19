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
#SBATCH --output=ADBDR_RunTappas.o
#SBATCH --error=ADBDR_RunTappas.e


# 04/11/2021: Run Kallisto of the RNA-Seq on the Iso-Seq scaffold (using the --rf-stranded)



##################################################################################################
#************************************* Prepare files for TappAS
# 6 files:
# 1) Iso-Seq Annotation scaffold file                       = WholeIsoSeq.collapsed.gff3
# 2) Iso-Seq retained isoforms from tama filtering          = tama_retained_pbid.txt
# 3) Iso-Seq Expression File
# 4) Iso-Seq Phenotype File
# 5) RNA-Seq Expression File
# 6) RNA-Seq Phenotype File


# Targeted Transcriptome: RNA-Seq alignment to Whole + Targeted Transcriptome
# 8) RNA-Seq Expression file - Kallisto RNA-Seq alignment to tama merged transcriptome of whole and targeted (all isoforms)
# 9) Annotation file - SQANTI3 annotation from tama merged transcriptome
# 10) List of Isoforms from target genes

# File 1,2
cp $PostIsoseq3_WKD/SQANTI3/$dataset".collapsed.gff3" $DiffAnalysis_WKD/TAPPAS_INPUT
cp $PostIsoseq3_WKD/SQANTI3_TAMA_FILTER/tama_retained_pbid.txt $DiffAnalysis_WKD/TAPPAS_INPUT

# File 3, 4, 7 - using Iso-Seq Expression
# counts_subset_4tappas <input_class> <output_class> <type_genes>
# regenerate expression matrix for tappas from long reads FL read counts
counts_subset_4tappas $PostIsoseq3_WKD/SQANTI3_TAMA_FILTER/$dataset"_sqantitamafiltered.classification.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression/$dataset"_sqantitamafiltered.expression.txt" AD
cp $RAWDIR/ADBDR_PhenotypeTAPPAS.txt $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression

# File 5, 6 - using RNA-Seq Expression
# Rscript script.R <input.dir> <output.file> <type=Whole/Targeted/WholeTargeted> <targeted.class.files>
source activate sqanti2_py3; Rscript $TAPPASFUNC/TAPPAS_RNASEQ_Exp.R $DiffAnalysis_WKD/RNASeq_SQANTI3 $dataset"RNASeq_sqantitamafiltered.expression.txt" Targeted $PostIsoseq3_WKD/SQANTI_TAMA_FILTER/$dataset"_sqantitamafiltered.classification.txt"
cp $DiffAnalysis_WKD/RNASeq_SQANTI3/AllBDRTargetedRNASeq_sqantitamafiltered.expression.txt $DiffAnalysis_WKD/TAPPAS_INPUT/RNASeq_Expression
cp $RAWDIR/ADBDR_RNASeqPhenotypeTAPPAS.txt $DiffAnalysis_WKD/TAPPAS_INPUT/RNASeq_Expression

##################################################################################################
#************************************* Run TappAS on Knight
scp -r sl693@login.isca.ex.ac.uk:/lustre/projects/Research_Project-MRC148213/lsl693/IsoSeq/Targeted_Transcriptome/ADBDR/Differential/TAPPAS_INPUT/* /mnt/data1/Szi/TAPPAS_ADBDR
/mnt/data1/Aaron/sw/jre1.8.0_181/bin/java -jar /mnt/data1/Szi/tappAS.jar

#### Tappas output
### tappAS projects
# Project.02103392608.tappas = Targeted ADBDR + IsoSeq
# Project.1802471916.tappas = Targeted ADBDR + RNASeq

cd $DiffAnalysis_WKD; mkdir TAPPAS_OUTPUT
cd $DiffAnalysis_WKD/TAPPAS_OUTPUT; mkdir IsoSeq_Expression RNASeq_Expression 
scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/tappasWorkspace/Projects/Project.02103392608.tappas/* $DiffAnalysis_WKD/TAPPAS_OUTPUT/IsoSeq_Expression/
scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/tappasWorkspace/Projects/Project.01950582580.tappas/* $DiffAnalysis_WKD/TAPPAS_OUTPUT/RNASeq_Expression/

mkdir -p $TAPPAS_OUTPUT_DIR/D_ONT
scp -r sLeung@knight.ex.ac.uk:/mnt/data1/Szi/tappasWorkspace/Projects/Project.1738176673.tappas/* $TAPPAS_OUTPUT_DIR/D_ONT
 
  
  