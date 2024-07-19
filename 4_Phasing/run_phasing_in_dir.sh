#!/bin/bash
DIR=$1
PLOIDY=4
MIN_BQ=13

echo "I'm going to $DIR"
STRAND=`grep "ref_strand=" $DIR/config|awk -F'=' '{print $2}'`
echo "strand is $STRAND"
PBID=`grep ">" $DIR/fake.fasta|awk -F"_" '{print $2}'`
PBID="fake_$PBID"

cat <<EOM >$DIR/run.sh

module load Miniconda2
module load Mamba
export PYTHONPATH=":/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake:/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake:/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake"
phasingDir=/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake/phasing
humanReferenceFasta=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa

# 3. create mpileup

if compgen -G "*phased*" > /dev/null; then
    echo "Phasing complete!"
else
  
  source activate nanopore

  minimap2 -ax splice fake.fasta ccs.fastq > ccs.sam
  #samtools view -b ../../all.fake.shortread.sorted.bam "fake_PB.10001" > mapped.shortread.bam
  
  # 3. create mpileup
  samtools view -bS ccs.sam > ccs.bam
  samtools sort ccs.bam > ccs.sorted.bam
  samtools mpileup --min-BQ 13 -f fake.fasta -s ccs.sorted.bam > ccs.mpileup
  
  #samtools merge out.bam mapped.shortread.bam ccs.sorted.bam
  #samtools mpileup --min-BQ 13 -f fake.fasta -s out.bam > ccs.mpileup
  
  # 4. run phasing
  mamba activate Spatial
  python /lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake/phasing/run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --strand - -o phased.nopartial -n 4
  python /lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake/phasing/run_phaser.py ccs.fastq ccs.sam ccs.mpileup fake.read_stat.txt fake.mapping.txt --partial_ok --strand - -o phased.partial -n 4
  
fi
EOM
