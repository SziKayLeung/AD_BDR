dat <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals/asts.tsv_asts_quant.tsv", header = T)
dat <- dat %>% mutate(transcript = word(transcript,c(1), sep = fixed("|")))
write.table(dat, "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals/asts.tsv_asts_quant_mod.tsv", sep = "\t", quote = F, row.names = F)
