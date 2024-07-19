#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=100:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address

gencode_gtf=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/gencode.v40.annotation.gtf
gencode_fa=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/mm10.fa
gencode_transcript_fasta=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/gencode.v38.transcripts.fa
gencode_translation_fasta=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/gencode.v38.pc_translations.fa
hexamer=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CPAT/Human_Hexamer.tsv
logit_model=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/CPAT/Human_logitModel.RData
classification=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/SQANTI3/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt
sqanti_fasta=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/SQANTI3/AllBDRTargeted.collapsed_classification.filtered_lite.fasta
sqanti_gtf=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/SQANTI3/AllBDRTargeted.collapsed_classification.filtered_lite.gtf
metamorpheus_toml=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/Task1SearchTaskconfig_orf.toml
meta_rescue_toml=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/Task1SearchTaskconfig_rescue_resolve.toml
mass_spec=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics_Rawdata
uniprot_fasta=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/reference_2019/Homo_sapiens.GRCh38.pep.all.fa
WKD=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics
name=ADBDR
coding_score_cutoff=0.0
min_junctions_after_stop_codon=2
lower_kb=1
upper_kb=4
lower_cpm=3

module load Miniconda2 
source activate lrp
export PATH=$PATH:/lustre/home/sl693/.dotnetexport PATH=$PATH:/lustre/home/sl693/.dotnet

# run_metamorpheus <input_fasta>
run_metamorpheus(){
  # variables 
  input_fasta=$1
  protein_dir=$2
  output_dir=$3
  output_name=$4

  echo "Processing $input_fasta"
  meta_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/softwares/MetaMorpheus/
  protein_data=$(for i in $protein_dir/*raw*; do echo $i; done)
  echo $protein_data

  cd $output_dir 
  dotnet $meta_dir/CMD.dll -g -o ./toml --mmsettings ./settings
  dotnet $meta_dir/CMD.dll -d $input_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $protein_data -t $metamorpheus_toml -v normal --mmsettings settings -o ./$output_name"_search_results"
  mv $output_name"_search_results"/Task1SearchTask/AllPeptides.psmtsv $output_name"_search_results"/Task1SearchTask/AllPeptides.$output_name".psmtsv"
  mv $output_name"_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.tsv $output_name"_search_results"/Task1SearchTask/AllQuantifiedProteinGroups.$output_name".filtered.tsv"  
}

protein_wkd=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics 
protein_dir=/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/Proteomics_Rawdata
cd $protein_wkd; mkdir -p All 
#run_metamorpheus $protein_wkd/$name".filtered_protein.fasta" $protein_dir $protein_wkd/All $name"_filtered"
#run_metamorpheus $protein_wkd/$name".protein_refined.fasta" $protein_dir $protein_wkd/All $name"_refined"
#run_metamorpheus $protein_wkd/$name"_hybrid.fasta" $protein_dir $protein_wkd/All $name"_hybrid"
run_metamorpheus $WKD/gencode_protein.fasta $protein_dir $protein_wkd/All Gencode
run_metamorpheus $uniprot_fasta $protein_dir $protein_wkd/All UniProt
