#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=1:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --array=0-2 # 3 samples
#SBATCH --output=2_IsoSeqPhasing-%A_%a.o
#SBATCH --error=2_IsoSeqPhasings-%A_%a.e


softwareDir=/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake/phasing/utils
humanReferenceFasta=/lustre/projects/Research_Project-MRC148213/lsl693/references/human/hg38.fa
cupcakeDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu
refineDir=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/3_refine
scriptDir=/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR/4_Phasing

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake
export PYTHONPATH="${PYTHONPATH}:/lustre/projects/Research_Project-MRC148213/lsl693/software/cDNA_Cupcake"

#cp ${softwareDir}/run_phasing_in_dir.sh .

cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing
ls -1d by_loci/*size*|xargs -n1 -i echo "bash ${scriptDir}/run_phasing_in_dir.sh {}" > cmd
bash cmd

ls -1d by_loci/*size*|xargs -n1 -i echo "cd {}; bash run.sh; cd ../../" > cmd2
#bash cmd2
#split -l 500 cmd2 splitcmd

files=($(ls *split*))
file=${files[${SLURM_ARRAY_TASK_ID}]}
bash ${file}
