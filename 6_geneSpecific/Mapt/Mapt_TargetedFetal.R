library(matrixStats)

cd /lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/Human_Mapt

export PATH=$PATH:/lustre/projects/Research_Project-MRC148213/lsl693/softwares/gffcompare

Fetal_Targeted="/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/Targeted/P0052_20220421_10661/Pool1/20220421_1656_1F_PAI87431_1fc2cdd3/P0052_analysis/Targeted_SQANTI3/Targeted_SQANTI3_classification.filtered_lite.gtf"

BDR_Targeted="/lustre/projects/Research_Project-MRC148213/lsl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/SQANTI3/AllBDRTargeted.collapsed_classification.filtered_lite.gtf"

cp $Fetal_Targeted .
cp $BDR_Targeted . 


gffcompare -r AllBDRTargeted.collapsed_classification.filtered_lite.gtf Targeted_SQANTI3_classification.filtered_lite.gtf -o BDRFetalTargeted


#### R
tmap = read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/Human_Mapt/BDRFetalTargeted.Targeted_SQANTI3_classification.filtered_lite.gtf.tmap", header= T)

tmap_exact <- tmap[tmap$class_code == "=",]
tmap_exact_mapt <- tmap[tmap$qry_gene_id == "ENSG00000186868.17",]

BDR_class <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/IsoSeq/Targeted_Transcriptome/ADBDR/Post_IsoSeq/SQANTI3/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt", header = T)

BDR_class$associated_transcript <- as.character(BDR_class$associated_transcript)
BDR_class["associated_transcript"][BDR_class["isoform"] == "PB.4515.63"] <- "ENST00000535772.6"
BDR_class["structural_category"][BDR_class["isoform"] == "PB.4515.63"] <- "full-splice_match"

BDR_detected <- BDR_class[BDR_class$associated_gene == "MAPT",] %>% 
  select(isoform, structural_category, associated_gene, associated_transcript, starts_with("FL.")) %>% 
  full_join(., tmap_exact_mapt[,c("ref_id","class_code","qry_id")],by=c("isoform"="ref_id")) %>% 
  mutate(detected = ifelse(is.na(class_code),"Unique","Common")) %>% 
  mutate(Sum = rowSums(as.matrix(select(., starts_with("FL."))), na.rm = TRUE)) 

ggplot(BDR_detected, aes(x = detected,y=Sum,colour = structural_category)) + geom_boxplot() + 
  geom_point(aes(fill = structural_category), size = 2, shape = 21, position = position_jitterdodge()) +
  labs(x = "Detection Status", y = "Sum of FL reads - BDR Targeted Dataset") + 
  theme_classic()


BDR_detected[BDR_detected$associated_gene == "MAPT" & BDR_detected$structural_category == "full-splice_match" & BDR_detected$detected == "Unique",]
BDR_detected[BDR_detected$associated_gene == "MAPT" & BDR_detected$structural_category == "full-splice_match" & BDR_detected$detected == "Common",]


