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

#************************************* DEFINE GLOBAL VARIABLES
# setting names of directory outputs
DiffAnalysis_WKD=/lustre/projects/Research_Project-MRC148213/lsl693/IsoSeq/Targeted_Transcriptome/ADBDR/Differential
PostIsoseq3_WKD=/lustre/projects/Research_Project-MRC148213/lsl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq
RNASeq_WKD=/lustre/projects/Research_Project-MRC148213/lsl693/IsoSeq/Targeted_Transcriptome/ADBDR/RNASeq
REFERENCE=/lustre/projects/Research_Project-MRC148213/lsl693/reference_2019
RNASeq_Filtered=/lustre/projects/Research_Project-193356/Project_10202/11_fastp_trimmed
DIFF_FUNC=/lustre/projects/Research_Project-MRC148213/lsl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis
TAPPASFUNC=/lustre/projects/Research_Project-MRC148213/lsl693/Scripts/IsoSeq_Tg4510/3_Differential_Analysis
RAWDIR=/lustre/projects/Research_Project-MRC148213/lsl693/Scripts/AD_BDR/Raw_Data

#cd $DiffAnalysis_WKD; mkdir RNASeq_SQANTI3 TAPPAS_INPUT
cd $DiffAnalysis_WKD/TAPPAS_INPUT; mkdir RNASeq_Expression IsoSeq_Expression

SAMPLES_NAMES=(BBN006_30024 BBN_24548 BBN002_28350 BBN_20194 BBN_24260 BBN003_30195 BBN002_29416 BBN_24287 BBN002_30035 BBN_15250 BBN002_29087 BBN003_26927 BBN_24289 BBN006_29902 BBN_25890 BBN002_29471 BBN_15247 BBN_18399 BBN002_26311 BBN_18405 BBN_24938 BBN006_29162 BBN_24253 BBN_19616 BBN_4240 BBN002_30920 BBN_15237 BBN_19632 BBN004_26227 BBN006_26347)
SAMPLE=${SAMPLES_NAMES[${SLURM_ARRAY_TASK_ID}]}

# sourcing functions script and input directories
FUNCTIONS=/lustre/projects/Research_Project-MRC148213/lsl693/Scripts/AD_BDR/1_Transcriptome_Annotation
source $FUNCTIONS/Targeted_ADBDR_Functions.sh

module load Miniconda2/4.3.21

dataset=AllBDRTargeted
sqname=$dataset".collapsed_classification.filtered_lite"

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


##################################################################################################
echo "#************************************* Prepare files for TappAS with further collapsing"
# counts_subset_4tappas <input_class> <output_class> <type_genes>
# regenerate expression matrix for tappas from long reads FL read counts
#counts_subset_4tappas $PostIsoseq3_WKD/SQANTI3_TAMA_FILTER/$dataset"_sqantitamafiltered.classification.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/$dataset"_sqantitamafiltered.expression.txt" AD

# Script <sqanti_filered_inputfile> <tama_filtered_inputfile> <output_expression_file> <output_finalisoform_file> <type == "Targeted/Whole">
#Rscript $DIFF_FUNC/Sqanti_Collapsed_Counts.R $PostIsoseq3_WKD/SQANTI3/$sqname"_classification.txt" $PostIsoseq3_WKD/SQANTI3/$dataset"_sqantitamafiltered.classification.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression/$dataset"_sqantisubset.expression" $DiffAnalysis_WKD/TAPPAS_INPUT/Retained_collapsed_pbid Targeted

#Rscript $DIFF_FUNC/Sqanti_Collapsed_Counts.R $PostIsoseq3_WKD/SQANTI3/$sqname"_classification.txt" $PostIsoseq3_WKD/SQANTI3/$dataset"_sqantitamafiltered.classification.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/IsoSeq_Expression/$dataset"_allsqantisubset.expression.txt" $DiffAnalysis_WKD/TAPPAS_INPUT/tama_retained_allcollapsed_pbid.txt All

## 12) collapse_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
# rerun on command line though logged files
# note tama_retained_collapsed_pbid.txt same as Retained_collapsed_pbid_FSMbylength.txt
#collapse_filter $DiffAnalysis_WKD/TAPPAS_INPUT/tama_retained_collapsed_pbid.txt $PostIsoseq3_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" $dataset $DiffAnalysis_WKD/COLLAPSE_FILTER/Length
#collapse_filter $DiffAnalysis_WKD/TAPPAS_INPUT/Retained_collapsed_pbid_FSMbyexpression.txt $PostIsoseq3_WKD/SQANTI3 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt" $dataset $DiffAnalysis_WKD/COLLAPSE_FILTER/Expression
