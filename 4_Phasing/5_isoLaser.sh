# /lustre/home/sl693/.local/lib/python3.10/site-packages/isolaser
SQANTI=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite.fasta 
SQANTIGTF=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite.gtf
SQANTICLASSFILE=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt

export REFERENCE=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/references
export GENOME_FASTA=$REFERENCE/human/hg38.fa
export GENOME_GTF=$REFERENCE/annotation/gencode.v40.annotation.gtf
export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/miscellaneous
export PYTHONPATH="${PYTHONPATH}:/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/miscellaneous"

tabix -p ${SQANTI}
isolaser_convert_gtf_to_fasta -g ${SQANTI} -f ${GENOME_FASTA} -o ADBDR

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/software/isoLASER
export PYTHONPATH="${PYTHONPATH}:/lustre/projects/Research_Project-MRC148213/lsl693/software/isoLASER"

cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/isoLaser
TRANSCRIPTOME_FASTA=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/gencode.v40.transcripts.fa
minimap2 -t 16 -ax splice:hq -uf --MD ${TRANSCRIPTOME_FASTA} ${SQANTI} > AllBDRTargeted.collapsed_classification.filtered_lite.sam
samtools view -bS AllBDRTargeted.collapsed_classification.filtered_lite.sam > AllBDRTargeted.collapsed_classification.filtered_lite.bam

#isolaser_filter_and_annotate -b {input.bam} -t AllBDRTargeted.collapsed_classification.filtered_lite.sam -g ${SQANTIGTF} -o {input.annot}
isoLaserDir=/lustre/projects/Research_Project-MRC148213/lsl693/software/isoLASER/src/isolaser/
python ${isoLaserDir}/annotate_reads.py -b AllBDRTargeted.collapsed_classification.filtered_lite.bam -t AllBDRTargeted.collapsed_classification.filtered_lite.sam  -g ${SQANTIGTF} -o input.annot.bam
isolaser_extract_exon_parts -g ${SQANTIGTF} -o transcript.db

source activate nanopore
cp ${SQANTIGTF} .
include_geneFeature_gtf.py -g AllBDRTargeted.collapsed_classification.filtered_lite.gtf -c ${SQANTICLASSFILE}

sed -i 's/transcript_id "NaN";//g' AllBDRTargeted.collapsed_classification.filtered_lite_geneIncluded.gtf 
sort -k1,1 -k4,4n AllBDRTargeted.collapsed_classification.filtered_lite_geneIncluded.gtf > AllBDRTargeted.collapsed_classification.filtered_lite_geneIncluded_sorted.gtf
bgzip  AllBDRTargeted.collapsed_classification.filtered_lite_geneIncluded_sorted.gtf
tabix -p gff AllBDRTargeted.collapsed_classification.filtered_lite_geneIncluded_sorted.gtf.gz

mamba activate Spatial
mkdir transcript.db
isolaser_extract_exon_parts -g AllBDRTargeted.collapsed_classification.filtered_lite_geneIncluded_sorted.gtf.gz -o transcript.db
