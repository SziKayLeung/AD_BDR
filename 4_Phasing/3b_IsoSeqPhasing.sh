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
#SBATCH --output=3b_IsoSeqPhasing.o
#SBATCH --error=3b_IsoSeqPhasings.e

module load Miniconda2 

softwareDir=/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake/phasing/utils
humanReferenceFasta=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
cupcakeDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu
refineDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/3_refine
scriptDir=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR/4_Phasing

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake
export PYTHONPATH="${PYTHONPATH}:/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake"

cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing

echo "Tagging*****************************"
awk -F "\t" '{ print $1 }' phased.no.partialID.txt  | tail -n +2 | xargs -n 1 -I {} bash -c '
tagPhasing(){
  echo "Tagging $1"
  cd $1
  filetype=(phased.no.partial phased.partial)
  softwareDir=/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake/phasing/utils
  
  for i in ${filetype[@]}; do 
    echo $i
    if [ -f ${i}.cleaned.hap_info.txt ]; then
      source activate sqanti2_py3
      python ${softwareDir}/tag_bam_post_phasing.py ccs.sorted.hg38.bam ${i}.cleaned.hap_info.txt ccs.sorted.hg38.${i}.tagged.bam
      
      source activate nanopore
      samtools index ccs.sorted.hg38.${i}.tagged.bam ccs.sorted.hg38.${i}.tagged.bam.bai
      samtools view -h -o ccs.sorted.hg38.${i}.tagged.sam ccs.sorted.hg38.${i}.tagged.bam
   else
      echo "${i}.cleaned.hap_info.txt not present"
    fi
  done 
  
}
tagPhasing {}
'
