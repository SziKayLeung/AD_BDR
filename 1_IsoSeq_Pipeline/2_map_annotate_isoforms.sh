#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=Targeted_ADBDR_Part2c.o
#SBATCH --error=Targeted_ADBDR_Part2c.e

# 11/05/2021: Run pipeline from refine to sqanti and tama (output: Targeted_ADBDR_Part2.o)
# 03/11/2021: Run kallisto with ADBDR RNA-Seq (30 samples), SQANTI3 with RNA-Seq support and TAMA filter (output: Targeted_ADBDR_Part2b.o)

##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
source $SC_ROOT/bdr_isoseq_config.sh
source $SC_ROOT/01_source_functions.sh


##-------------------------------------------------------------------------
echo "#*************************************  Isoseq3 [Function 1, 2, 3, 6]"
echo "Already processed in batch (Mouse/Targeted_Mouse_Part1.sh) Functions 1 - 3, 6 for all batched runs"

## 4) run_targeted_REFINE 
run_targeted_REFINE 

# 6) run_CLUSTER $Sample $Input_REFINE_directory $Output_directory
for i in ${ALL_SAMPLES_NAMES[@]}; do run_CLUSTER $i; done

##-------------------------------------------------------------------------
echo "#************************************* Isoseq3 and Post_Isoseq3 [Function 5, 7, 8]"
## 5) merging_at_refine <output_name> <samples.....>
merging_at_refine $NAME ${ALL_SAMPLES_NAMES[@]}

## 7) run_map_cupcakecollapse <output_name> 
run_map_cupcakecollapse $NAME 

## 8) demux_targeted <refine_dir> <input_cluster_report> <input_tofu_readstat> <output_path_file>
demux_targeted $NAME

##-------------------------------------------------------------------------
echo "#************************************* RNAseq [Function 9, 10]"
## 9) run_star <list_of_samples> <input_directory> <output_dir>
bash $SC_ROOT/2b_map_rnaseq_genome.sh
wait
# .e files were empty; cat all output files > ADBDR_Rnaseq_Star.o

## 10) rnaseq_merge_fastq <sample>
rnaseq_merge_fastq $NAME

##-------------------------------------------------------------------------
echo "#************************************* RNASeq & IsoSeq [Function 11]"
## 11) run_kallisto <sample>
run_kallisto $NAME 

##-------------------------------------------------------------------------
echo "#************************************* SQANTI3 [Function 12]"
## 12) run_sqanti3 <sample> <mode=basic/full/nokallisto/lncrna>
run_sqanti3 $NAME full
run_sqanti3 $NAME basic

##-------------------------------------------------------------------------
echo "#************************************* TAMA filter [Function 13,14,16]"
## 13) TAMA_remove_fragments <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_remove_fragments $NAME full

## 14) TAMA_sqanti_filter <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_sqanti_filter $NAME full

## 16) TAMA_tappas_input <sample> <mode=basic/full/nokallisto/lncrna>
TAMA_tappas_input $NAME full

##-------------------------------------------------------------------------
echo "#************************************* QC [Function 15]"
# parse_stats_per_sample <sample>
parse_stats_per_sample $NAME

run_target_rate


##-------------------------------------------------------------------------
echo "#************************************* Characterisation of isoforms"

# run_cpat <input_fasta> <output_name>
run_cpat $WKD_ROOT/9_sqanti3/basic/$NAME.collapsed_corrected.fasta $NAME

# extract_best_orf <sample> 
extract_best_orf $NAME

# remove_3ISM <sample> <mode=basic/full/nokallisto/lncrna>
remove_3ISM $NAME basic

# colour_by_abundance <sample> <input_gtf>
colour_by_abundance $NAME $WKD_ROOT/9_sqanti3/basic/$NAME.collapsed_corrected.gtf 

# subset_gene_reference 
subset_gene_reference 
#extract_reference_info.py

# full_characterisation <gene> 
full_characterisation $NAME SNCA

