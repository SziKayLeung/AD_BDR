
module load Mamba
mamba activate Spatial

extractScript=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/phasing/extract_variant_per_gene.py 
dir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/genotype
phenotype=BDR_pathology_data_for_analysis.csv

TGENES=(ABCA1 SORL1 MAPT BIN1 TARDBP APP ABCA7 PTK2B ANK1 FYN CLU CD33 FUS PICALM SNCA APOE TRPA1 RHBDF2 TREM2 VGF)
for gene in ${TGENES[@]}; do 
  echo ${gene}
  python ${extractScript} ${dir}/NeuroChipDiseaseVariants.txt --gene ${gene} -b BDR_imputed_EUR_QCd --dir ${dir}/imputed -p ${dir}/BDR_pathology_data_for_analysis.csv
done