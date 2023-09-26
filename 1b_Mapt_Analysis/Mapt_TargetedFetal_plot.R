# identify human tau isoforms in mouse targeted transcriptome dataset 
human_dir = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/ONT/Targeted_Transcriptome/TALON_Human/Filtered"
human.class.files = SQANTI_class_preparation(paste0(human_dir, "/SQANTI3/ONTTargeted_filtered_talon_classification.filtered_lite_classification.txt"),"nstandard") %>% filter(associated_gene == "MAPT")
human.abundance = read.table(paste0(human_dir, "/ONTTargeted_filtered_talon_abundance_filtered.tsv"), sep = "\t", header = T)
human.class.files = merge(human.class.files, human.abundance, by.x = "isoform", by.y = "annot_transcript_id") 
targeted.phenotype = read.table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/IsoSeq_Tg4510/Raw_Data/Targeted_Transcriptome/TargetedMouse_PhenotypeTAPPAS.txt", header = T) %>% 
  mutate(sampleID = paste0(word(sample, c(2), sep = fixed("."))),
         genotype = ifelse(group == "CONTROL", "WT","TG"),
         age = paste0(time,"mos")) %>% 
  mutate(colID = paste0(sampleID,"_",genotype,"_",age))

samples = c("S19","K24","L22","M21","O18","O23","O22","P19","T20","Q20","Q21","S18","S23","Q18","Q17","L18","Q23","T18")                  
df = human.class.files %>% select(isoform, associated_transcript, associated_gene, samples) %>% reshape2::melt() %>% 
  left_join(., targeted.phenotype, by=c("variable"="sampleID")) 

plot_ont_expression <- function(isoform){
  p <- ggplot(df[df$isoform == isoform,], aes(x = as.factor(time), y = value, colour = group)) + 
    geom_boxplot() + geom_point(position=position_jitterdodge()) + mytheme +
    labs(x = "Age (months)", y = "ONT FL reads", title = isoform) + 
    scale_colour_manual(values = c(label_colour("TG"),label_colour("WT")), labels=c("TG","WT")) + 
    theme(legend.position = "bottom")
  
  return(p)
}

plot_ont_expression("TALONT000279659")
plot_ont_expression("TALONT000280246")
plot_ont_expression("ENST00000351559.10")
plot_ont_expression("ENST00000334239.12")


ggplot(df, aes(x = as.factor(time), y = value, colour = isoform)) + geom_point() + facet_grid(~genotype) + 
  labs(x = "Age (months)", y = "Expression")

df[df$isoform == "TALONT000279659",]
find_mapt_ont("human")
