#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=12 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=4a_run_talon_database.o
#SBATCH --error=4a_run_talon_database.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config
source $SC_ROOT/1_ONT_Pipeline/01_source_functions.sh


##-------------------------------------------------------------------------

if [ ! -f ${TALON_DB}/${SPECIES}_talon.db ]; then
  
  source activate sqanti2_py3
  
  echo "Generating database: ${SPECIES}_talon"
  cd ${TALON_DB}
  talon_initialize_database --f ${GENOME_GTF} --a ${SPECIES}_annot --g ${SPECIES} --o ${SPECIES}_talon &> ${SPECIES}_talon_init.log
  
  source deactivate

else
  echo "Database present: ${SPECIES}_talon"

fi
