module load Miniconda2/4.3.21
source activate ficle
FICLE_ROOT=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/FICLE/
export PATH=$PATH:${FICLE_ROOT}

dir=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/6_sqanti/sqanti
output_dir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/devFicle
cpat=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARI/7_cpat/cpat
GENOME_GTF=/lustre/projects/Research_Project-MRC148213/lsl693/references/annotation/gencode.v40.annotation.20Targets.gtf
grep ONT17_1068 ${dir}/WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf > ${output_dir}/MAPT_WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf


ficle.py --genename=MAPT --geneid=ONT17_1068 \
--reference=${GENOME_GTF} \
--input_gtf=${output_dir}/MAPT_WholeTargeted_cleaned_aligned_merged_collapsed_qced_corrected_2reads2samples_2reads2samples_nomonointergenic.gtf \
--input_class=${output_dir}/Mapt_WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_2reads2samples_classification_noMonoIntergenic.txt \
--cpat=${cpat}/WholeTargeted_fixed.ORF_prob.best.tsv   \
--output_dir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/devFicle 
