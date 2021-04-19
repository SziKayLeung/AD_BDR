#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=10:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Targeted_ADBDR_Part2.o
#SBATCH --error=Targeted_ADBDR_Part2.e

#************************************* DEFINE GLOBAL VARIABLES
# File directories
Isoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/IsoSeq
PostIsoseq3_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq
RNASeq_WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/RNASeq

# For Pooled Targeted
### Important order of BAM files is the same as the sample names
BATCH_NAMES=(Targeted_Seq_1)
RAWDIR=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/1_Transcriptome_Annotation/Raw_Data

# Other input files and directory
REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019

# For Demultiplexing Samples
Pooled_Samples_Targeted1=(A1 A2 A3 A4 A5 A6 A7 A8 A9 A10)
Barcoded_Targeted1_config_file=$RAWDIR/Barcode_Configs/Isoseq_ADBDR_Targeted1_barcode.config
ALL_SAMPLES_NAMES=(A1 A2 A3 A4 A5 A6 A7 A8 A9 A10)

# sourcing functions script
FUNCTIONS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/1_Transcriptome_Annotation
source $FUNCTIONS/Targeted_ADBDR_Functions.sh

module load Miniconda2/4.3.21
################################################################################################
echo "#*************************************  Isoseq3 [Function 1, 2, 3, 6]"
echo "Already processed in batch (Mouse/Targeted_Mouse_Part1.sh) Functions 1 - 3, 6 for all batched runs"

## 4) run_targeted_REFINE ${BATCH} $Input_config_file $Input_Lima_sample $Input_LIMA_directory $Output_directory
for i in ${Pooled_Samples_Targeted1[@]}; do run_targeted_REFINE $i $Barcoded_Targeted1_config_file Targeted_Seq_1 $Isoseq3_WKD/LIMA $Isoseq3_WKD/REFINE; done

# 6) run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
for i in ${ALL_SAMPLES_NAMES[@]}; do run_CLUSTER $i $Isoseq3_WKD/REFINE $Isoseq3_WKD/CLUSTER; done

################################################################################################
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8]"
## 5) merging_at_refine <input_flnc_bam_dir> <output_directory> <output_name> <samples.....>
# Targeted_Seq_3b = complete run of Batch 3
merging_at_refine $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER AllBDRTargeted Targeted_Seq_1

## 7) run_map_cupcakecollapse <sample_prefix_input/output_name> <isoseq3_input_directory> <mapping_output_directory> <tofu_output_directory>
run_map_cupcakecollapse AllBDRTargeted $Isoseq3_WKD/MERGED_CLUSTER $PostIsoseq3_WKD/MAP $PostIsoseq3_WKD/TOFU

## 8) demux_targeted <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
demux $Isoseq3_WKD/REFINE $Isoseq3_WKD/MERGED_CLUSTER/AllBDRTargeted.clustered.cluster_report.csv $PostIsoseq3_WKD/TOFU/AllBDRTargeted.collapsed.read_stat.txt $PostIsoseq3_WKD/TOFU/AllBDRTargeted.Demultiplexed_Abundance.txt

################################################################################################
echo "#************************************* RNAseq [Function 9, 10]"
## 9) run_star <list_of_samples> <J20/Tg4510_input_directory> <output_dir>

## 10) mouse_merge_fastq <RNASEQ_input_dir> <Kallisto_output_dir> <sample_prefix_output_name>

################################################################################################
echo "#************************************* RNASeq & IsoSeq [Function 11]"
## 11) run_kallisto <sample_prefix_output_name> <input_tofu_fasta> <merged_fastq_input_dir> <output_dir>

################################################################################################
echo "#************************************* SQANTI2 [Function 12]"
## 12) run_sqanti2 <input_tofu_prefix> <input_gtf> <input_tofu_dir> <input_RNASEQ_dir> <input_KALLISTO_file> <input_abundance> <output_dir> <mode=genome/noexp/lncrna>
run_sqanti2 AllBDRTargeted.collapsed AllBDRTargeted.collapsed.gff $PostIsoseq3_WKD/TOFU $RNASeq_WKD/MAPPED $PostIsoseq3_WKD/KALLISTO/AllBDRTargeted.mod.abundance.tsv $PostIsoseq3_WKD/TOFU/AllBDRTargeted.Demultiplexed_Abundance.txt $PostIsoseq3_WKD/SQANTI2 genome

################################################################################################
echo "#************************************* TAMA filter [Function 13,14]"
## 13) TAMA_remove_fragments <input_collapsed.filtered.gff> <input/output_prefix_name> <input/output_dir>
TAMA_remove_fragments $PostIsoseq3_WKD/SQANTI2/AllBDRTargeted.collapsed_classification.filtered_lite.gtfAllBDRTargeted $PostIsoseq3_WKD/TAMA

## 14) TAMA_sqanti_filter <TAMA_remove_fragments.output> <sqanti_filtered_dir> <sqanti_output_txt> <sqanti_output_gtf> <sqanti_output_fasta> <output_prefix_name> <output_dir>
sqname=AllBDRTargeted.collapsed_classification.filtered_lite
TAMA_sqanti_filter $PostIsoseq3_WKD/TAMA/AllBDRTargeted.bed $PostIsoseq3_WKD/SQANTI2 $sqname"_classification.txt" $sqname".gtf" $sqname".fasta" $sqname"_junctions.txt"AllBDRTargeted $PostIsoseq3_WKD/SQANTI_TAMA_FILTER

################################################################################################
echo "#************************************* QC [Function 15]"
# parse_stats_per_sample <input_ccs.bam_dir> <Input_LIMA_directory> <output_prefix_name>
parse_stats_per_sample $Isoseq3_WKD/CCS $Isoseq3_WKD/LIMA AllBDRTargeted
