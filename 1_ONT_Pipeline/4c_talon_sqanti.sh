#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mem=200G # specify bytes of memory to reserve
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=4c_talon_sqanti.o
#SBATCH --error=4c_talon_sqanti.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------


# generate config file 
#ls ${WKD_ROOT}/5_talon_label/*labeled.sam* > ${WKD_ROOT}/6_talon/labelled_files.csv
#source activate sqanti2_py3
#Rscript ${TALONCONFIGGEN} ${BARCODE_CONFIG} ${WKD_ROOT}/6_talon/labelled_files.csv 
#cat ${WKD_ROOT}/6_talon/talon.config

# run_talon <config_file> <output_dir>
#run_talon ${WKD_ROOT}/6_talon/talon.config
#post_talon ${WKD_ROOT}/6_talon

# convert talon expression file 
source activate nanopore 
Rscript ${TALEXP} -e ${UNFIL_DIR}/${NAME}_unfiltered_talon_abundance.tsv -o ${UNFIL_DIR} -n ${NAME}_unfiltered_talon_abundance_sqinput
Rscript ${TALEXP} -e ${FIL_DIR}/${NAME}_filtered_talon_abundance.tsv -o=${FIL_DIR} -n=${NAME}_filtered_talon_abundance_sqinput

# run sqanti3 
run_sqanti3 ${UNFIL_DIR}/${NAME}_unfiltered_talon.gtf ${WKD_ROOT}/7_sqanti3/1_unfiltered
run_sqanti3 ${FIL_DIR}/${NAME}_filtered_talon.gtf ${WKD_ROOT}/7_sqanti3/2_filtered

sqpath=${WKD_ROOT}/7_sqanti3/2_filtered/${NAME}_filtered_talon
source activate sqanti2_py3
cd ${WKD_ROOT}/7_sqanti3/2_filtered
python ${ISOANNOT} ${sqpath}_corrected.gtf ${sqpath}_classification.txt ${sqpath}_junctions.txt -gff3 ${GFF3}
