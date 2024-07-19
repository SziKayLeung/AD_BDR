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
#SBATCH --output=1_mass_spec_transfer.o
#SBATCH --error=1_mass_spec_transfer.e

# 01/08/2022: Transferred raw mass-spectrometry data for remaining 69 core AD-BDR samples (from filesender)

cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/Proteomics/1b_raw_remaining

wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690022" -O 20210108_AD_2.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690024" -O 20210108_AD_13.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690028" -O 20210108_AD_26.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690030" -O 20210108_AD_37.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690032" -O 20210108_AD_49.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690034" -O 20210108_AD_61.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690036" -O 20210108_AD_62.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690038" -O 20210108_AD_73.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690040" -O 20210108_AD_74.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690042" -O 20210108_AD_85.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690044" -O 20210108_AD_86.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690046" -O 20210110_AD_25.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690048" -O 20210110_AD_38.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690050" -O 20210110_AD_50.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690056" -O 20210112_AD_15.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690060" -O 20210112_AD_27.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690062" -O 20210112_AD_28.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690064" -O 20210112_AD_39.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690066" -O 20210112_AD_40.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690068" -O 20210112_AD_51.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690070" -O 20210112_AD_52.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690072" -O 20210112_AD_64.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690074" -O 20210112_AD_75.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690076" -O 20210112_AD_76.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690092" -O 20210115_AD_29.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690094" -O 20210115_AD_30.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690096" -O 20210115_AD_41.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690098" -O 20210115_AD_42.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690102" -O 20210115_AD_54.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690104" -O 20210115_AD_65.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690106" -O 20210115_AD_66.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690108" -O 20210115_AD_77.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690110" -O 20210115_AD_78.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690116" -O 20210118_AD_63.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690120" -O 20210120_AD_19.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690122" -O 20210120_AD_31.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690124" -O 20210120_AD_43.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690126" -O 20210120_AD_53.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690128" -O 20210120_AD_55.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690130" -O 20210120_AD_67.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690132" -O 20210120_AD_79.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690138" -O 20210308_AD_24.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690140" -O 20210308_AD_36.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690142" -O 20210308_AD_48.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690144" -O 20210308_AD_60.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690146" -O 20210308_AD_72.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690152" -O 20210324_AD_11.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690154" -O 20210324_AD_23.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690156" -O 20210324_AD_35.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690158" -O 20210324_AD_47.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690160" -O 20210324_AD_59.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690162" -O 20210324_AD_71.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690168" -O 20210325_AD_10.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690170" -O 20210325_AD_22.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690172" -O 20210325_AD_34.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690174" -O 20210325_AD_46.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690176" -O 20210325_AD_58.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690178" -O 20210325_AD_70.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690186" -O 20210328_AD_21.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690188" -O 20210328_AD_33.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690190" -O 20210328_AD_45.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690192" -O 20210328_AD_57.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690194" -O 20210328_AD_69.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690196" -O 20210328_AD_81.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690202" -O 20210330_AD_20.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690204" -O 20210330_AD_32.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690206" -O 20210330_AD_44.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690208" -O 20210330_AD_56.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690210" -O 20210330_AD_68.raw
wget "https://filesender.surf.nl/download.php?token=f391b59f-47e8-45de-b742-db6d99fcf6ef&files_ids=8690214" -O 20210330_AD_92.raw
