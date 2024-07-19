#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=S.K.Leung@exeter.ac.uk # email address
#SBATCH --output=4_runBambu.o
#SBATCH --error=4_runBambu.e


##-------------------------------------------------------------------------

# source config file and function script
module load Miniconda2/4.3.21
SC_ROOT=/lustre/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR
LOGEN_ROOT=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/LOGen
source $SC_ROOT/1_ONT_Pipeline/bdr_ont.config

# input variables 
alignedDir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/D_ONT/3_minimap/Batch2
annoRDS=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/0_metadata/Bambu_gencode_v40_annotations.rds
outputDir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/D_ONT/5_bambu
  
##-------------------------------------------------------------------------

source deactivate
module load R
Rscript ${LOGEN_ROOT}/assist_ont_processing/run_Bambu.R -i ${alignedDir} -f ${GENOME_FASTA} -a ${annoRDS} -o ${outputDir}


#test.bam <- list.files(path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/A_IsoSeq/6_minimap/Individual/", pattern = "bam", full = T)
#print(test.bam)
#test.bam <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/A_IsoSeq/6_minimap/Individual/A02.clustered.hq.fasta.bam"
#fa.file <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/references/human/hg38.fa"
#gtf.file <- "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/references/annotation/gencode.v40.annotation.gtf"

#annotations <- prepareAnnotations(gtf.file)
#saveRDS(annotations, "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/0_metadata/Bambu_gencode_v40_annotations.rds" )
#annotations <- readRDS("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/0_metadata/Bambu_gencode_v40_annotations.rds")
#se <- bambu(reads <- test.bam, annotations = annotations, genome = fa.file)
#writeBambuOutput(se, path = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/A_IsoSeq/7_bambu/")