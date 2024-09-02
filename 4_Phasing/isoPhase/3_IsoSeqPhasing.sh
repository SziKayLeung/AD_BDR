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
#SBATCH --output=3_IsoSeqPhasing.o
#SBATCH --error=3_IsoSeqPhasings.e

module load Miniconda2 

softwareDir=/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake/phasing/utils
humanReferenceFasta=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
cupcakeDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu
refineDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/3_refine
scriptDir=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR/4_Phasing

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake
export PYTHONPATH="${PYTHONPATH}:/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake"

cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing
##python ${softwareDir}/summarize_byloci_results.py
##python ${softwareDir}/collect_all_vcf.py --vcf phased.partial.cleaned.vcf


#awk -F "\t" '{ print $1 }' summarized.isophase_results.txt | tail -n +2 | xargs -n 1 -I {} bash -c '
#if [ -f {}/phased.nopartial.cleaned.hap_info.txt ]; then
#  echo {} >> phased.no.partialID.txt
#fi
#'

echo "Aligning*****************************"

awk -F "\t" '{ print $1 }' phased.no.partialID.txt | tail -n +2 | xargs -n 1 -I {} bash -c '
alignGenome(){
  humanReferenceFasta=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
  echo "Aligning $1 to genome"
  cd $1
  ## 5. create genome CCS BAM file for reference
  if [ ! -f ccs.sorted.hg38.bam ]; then
    source activate nanopore
    minimap2 -ax splice ${humanReferenceFasta} ccs.fasta > ccs.hg38.sam
    samtools view -bS ccs.hg38.sam > ccs.hg38.bam
    samtools sort ccs.hg38.bam > ccs.sorted.hg38.bam
    samtools index ccs.sorted.hg38.bam ccs.sorted.hg38.bam.bai
    echo "Alignment complete!"
  fi
}
alignGenome {}
'

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

mkdir finalbam
awk -F "\t" '{ print $1 }' summarized.isophase_results.txt | tail -n +2 | xargs -n 1 -I {} bash -c '
echo {}
sample=$(basename {})
if [ -f {}/ccs.sorted.hg38.phased.partial.tagged.bam ]; then
  cp {}/ccs.sorted.hg38.phased.partial.tagged.bam finalbam/${sample}ccs.sorted.hg38.phased.partial.tagged.bam
  cp {}/ccs.sorted.hg38.phased.partial.tagged.bam.bai finalbam/${sample}ccs.sorted.hg38.phased.partial.tagged.bam.bai
fi

if [ -f {}/ccs.sorted.hg38.phased.no.partial.tagged.bam ]; then
  cp {}/ccs.sorted.hg38.phased.no.partial.tagged.bam finalbam/${sample}ccs.sorted.hg38.phased.no.partial.tagged.bam
  cp {}/ccs.sorted.hg38.phased.no.partial.tagged.bam.bai finalbam/${sample}ccs.sorted.hg38.phased.no.partial.tagged.bam.bai
fi
'

phasedPartialTaggedBam=$(ls finalbam/*phased.partial.tagged.bam*)
source activate nanopore
samtools merge merged.phased.partial.tagged.bam ${phasedPartialTaggedBam}
samtools index merged.phased.partial.tagged.bam merged.phased.partial.tagged.bam.bai
# no phased.no.partial

# remove commands
rm *cmd*

cd finalbam
ls *bai* > ../finalTagged.txt
sed -i -e 's/ccs.sorted.hg38.phased.partial.tagged.bam.bai//g' ../finalTagged.txt

samtools view merged.phased.partial.tagged.bam > merged.phased.partial.tagged.sam