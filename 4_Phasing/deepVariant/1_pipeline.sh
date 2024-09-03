
# Download Precompiled Binaries
wget https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip
export PATH=/lustre/projects/Research_Project-MRC148213/lsl693/software/gatk-4.6.0.0/gatk:$PATH
export PATH=/lustre/projects/Research_Project-MRC148213/lsl693/software/java/jdk-21.0.4+7-jre/bin:$PATH

GENCODE_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
gatk CreateSequenceDictionary -R ${GENCODE_FASTA}

alignedDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/6_minimap
outputDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/deepVariant

gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=16" SplitNCigarReads \
-R ${GENCODE_FASTA} \
-I ${alignedDir}/AllBDRTargeted.sorted.bam \
-O ${outputDir}/AllBDRTargeted.sorted_sncr.bam


/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.dict