#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=S.K.Leung@exeter.ac.uk # email address
#SBATCH --output=convertFast52Pod5.o
#SBATCH --error=convertFast52Pod5.e

# raw data ONT batch2 to RDS 
#cp -r /lustre/recovered/Research_Project-MRC148213/sl693/AD_BDR/1_raw/* /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/1_raw/
module load Miniconda2
source activate sqanti2_py3

O100859=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/1_raw/B_ONT/P0063_20221026_100859/BDR_enriched/20221026_1304_2F_PAM95037_d0bc299c
cd ${O100859}
pod5 convert fast5 ./fast5_pass/*.fast5 --output pod5/ --one-to-one ./fast5_pass/
pod5 convert fast5 ./fast5_fail/*.fast5 --output pod5/ --one-to-one ./fast5_fail/

O10859=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/1_raw/B_ONT/P0063_20221031_10859/BDR_enriched/20221031_1830_2F_PAM95037_04a8d91e/
cd ${O10859}
pod5 convert fast5 ./fast5_pass/*.fast5 --output pod5/ --one-to-one ./fast5_pass/
pod5 convert fast5 ./fast5_fail/*.fast5 --output pod5/ --one-to-one ./fast5_fail/
  
O10918=/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/1_raw/B_ONT/P0073_20230209_10918/BDR/20230209_1549_3F_PAK63174_1ae1cd51/
cd ${O10918}
pod5 convert fast5 ./fast5_pass/*.fast5 --output pod5/ --one-to-one ./fast5_pass/
pod5 convert fast5 ./fast5_fail/*.fast5 --output pod5/ --one-to-one ./fast5_fail/
  
