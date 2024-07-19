protein_tsv = read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/Proteomics/ADBDR.sqanti_protein_classification.tsv", sep = "\t", header = T)

protein_tsv = merge(targeted.class.files[,c("isoform","associated_gene")], protein_tsv,by.x = "isoform",by.y="pb")

protein_tsv %>% group_by(associated_gene, pr_splice_cat) %>% tally() %>%
  ggplot(., aes(x = associated_gene, y = n, fill = pr_splice_cat)) + geom_bar(stat = "identity") + 
  labs(x = "Target Gene", y = "Number of protein isoforms") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_discrete(name = "Category") + mytheme

nrow(targeted.class.files[targeted.class.files$associated_gene == "BIN1",])
nrow(targeted.class.files[targeted.class.files$associated_gene == "BIN1" & targeted.class.files$coding == "coding",])
nrow(protein_tsv[protein_tsv$associated_gene == "BIN1",])

TargetGene <- c("ABCA1","SORL1","MAPT","BIN1","TARDBP","APP","ABCA7","PTK2B","ANK1","FYN","CLU","CD33","FUS","PICALM","SNCA","APOE","TRPA1","RHBDF2","TREM2","VGF")

output <- data.frame()
for(i in 1:length(TargetGene)){
  output[i, 1] <- TargetGene[i]
  output[i, 2] <- nrow(targeted.class.files[targeted.class.files$associated_gene == TargetGene[i],])
  output[i, 3] <- nrow(targeted.class.files[targeted.class.files$associated_gene == TargetGene[i] & targeted.class.files$coding == "coding",])
  output[i, 4] <- nrow(protein_tsv[protein_tsv$associated_gene == TargetGene[i],])
  output[i, 5] <- nrow(protein_tsv[protein_tsv$associated_gene == TargetGene[i] & protein_tsv$pr_splice_cat == "full-splice_match",])
  output[i, 6] <- nrow(protein_tsv[protein_tsv$associated_gene == TargetGene[i] & protein_tsv$pr_splice_cat == "incomplete-splice_match",])
  output[i, 7] <- nrow(protein_tsv[protein_tsv$associated_gene == TargetGene[i] & protein_tsv$pr_splice_cat == "novel_not_in_catalog",])
  output[i, 8] <- nrow(protein_tsv[protein_tsv$associated_gene == TargetGene[i] & protein_tsv$pr_splice_cat == "novel_in_catalog",])
}
colnames(output) <- c("associated_gene","total_num_transcripts","total_num_coding_transcripts","total_protein","pFSM","pISM","pNNC","pNIC")


reshape::melt(output) %>% filter(variable %in% c("total_num_transcripts","total_num_coding_transcripts","total_protein")) %>% 
  ggplot(., aes(x = associated_gene, y = value, fill = variable)) + geom_bar(stat = "identity", position = position_dodge()) +
  mytheme + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") + 
  labs(x = "Target Gene", y = "Number") + 
  scale_fill_discrete(name = "Category", labels = c("Transcripts","Coding Transcripts","Protein")) 

peptides <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/Proteomics/ADBDR_refined_peptides.bed12", header = F, sep = "\t")
peptide_tally <- peptides %>% mutate(gene = word(word(V4, c(2), sep = fixed("(")), c(1), sep = fixed(")"))) %>% 
  group_by(gene) %>% tally()
View(peptide_tally)

peptide_tally %>% filter(gene %in% TargetGene)
