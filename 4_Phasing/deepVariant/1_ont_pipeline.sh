#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=50:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.k.leung@exeter.ac.uk # email address
#SBATCH --output=deepVariant.o
#SBATCH --error=deepVariant.e

# run deepVariant on ONT BDR dataset
##-------------------------------------------------------------------------

export PATH=/lustre/projects/Research_Project-MRC148213/lsl693/software/udocker-1.3.17/udocker:$PATH
export PATH=/lustre/projects/Research_Project-MRC148213/lsl693/software/gatk-4.6.0.0:$PATH
export PATH=/lustre/projects/Research_Project-MRC148213/lsl693/software/java/jdk-21.0.4+7-jre/bin:$PATH
export SOFTWARE=/lustre/projects/Research_Project-MRC148213/lsl693/software
export THREADS=16
GENCODE_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
alignedDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/D_ONT/3_minimap/Batch2
outputDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/deepVariant


# ID
export METADATA=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/0_metadata/B_ONT/Selected_ONTTargeted_BDR.csv
ID=($(awk -F, '{OFS=",";print $2}' ${METADATA}))
BARCODE=($(awk -F, '{OFS=",";print $15}' ${METADATA}))
for i in {1..53}; do 
  #echo $i
  echo ${ID[$i]}
  echo Barcode: ${BARCODE[$i]}
  ls "${alignedDir}/BC${BARCODE[$i]}_merged_combined_minimap2.log"
done 

# RG tag: read group identifier, barcode
# SM tag: sample ID
# LB tag: library name
# PL tag: platform
samtools addreplacerg -r "@RG\tID:BC32\tSM:Sample1\tLB:lib1\tPL:ONT" -o output_with_rg.bam input.bam

##-------------------------------------------------------------------------
# Download Precompiled Binaries
#wget https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip
#gatk CreateSequenceDictionary -R ${GENCODE_FASTA}

cd ${outputDir}
#gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=16" SplitNCigarReads -R ${GENCODE_FASTA} -I ${alignedDir}/AllBDRTargeted.sorted.bam -O ${outputDir}/AllBDRTargeted.sorted_sncr.bam

#module load SAMtools/1.9-foss-2018b
#Rscript ${SOFTWARE}/lrRNAseqVariantCalling/tools/flagCorrection.r ${alignedDir}/AllBDRTargeted.sorted.bam ${outputDir}/AllBDRTargeted.sorted_sncr.bam ${outputDir}/AllBDRTargeted.sorted_sncr_fc.bam ${THREADS}

### index
samtools index -@ ${THREADS} ${outputDir}/AllBDRTargeted.sorted_sncr_fc.bam

# run deepVariant
INPUT_DIR=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/deepVariant
OUTPUT_DIR=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/deepVariant
REF_DIR=/lustre/projects/Research_Project-MRC148213/lsl693/references/human
BIN_VERSION="1.3.0"
udocker run \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${REF_DIR}:${REF_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=${REF_DIR}/hg38.fa \
  --reads=${INPUT_DIR}/AllBDRTargeted.sorted_sncr_fc.bam  \
  --output_vcf=${OUTPUT_DIR}/deepvariant_calls.vcf \
  --verbosity=3 --num_shards=${THREADS} --logging_dir=${OUTPUT_DIR}/logs 