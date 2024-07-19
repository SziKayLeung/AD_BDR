
suppressMessages(library("cowplot"))
redundant_iso = proteinInput$t.class.files[proteinInput$t.class.files$associated_gene == "APOE", "isoform"]
ApoeIsoList <- data.frame(
  Isoform = unlist(IsoList <- list(
    subset = redundant_iso,
    protein = as.character(unique(proteinInput$peptide_orf[proteinInput$peptide_orf$isoform %in% redundant_iso,"transcript_id"]))
  )),
  Category = rep(names(IsoList), lengths(IsoList))
)
ApoeIsoList$colour <- c(rep(NA,length(ApoeIsoList$Category[ApoeIsoList$Category != "DTE"])))

# Apoe
pApoeExp <- plot_trans_exp_individual("PB.5435.41",annoResTran$pWald$norm_counts %>% mutate(value = normalised_counts))
pApoeTracks <- ggTranPlots(proteinInput$merged_peptide_gtf, proteinInput$t.class.files, isoList = c(as.character(ApoeIsoList$Isoform)), selfDf = ApoeIsoList, gene = "APOE")
plot_grid(pApoeExp, pApoeTracks)

proteinInput$t.class.files[proteinInput$t.class.files$base_acc == "PB.5435.41",]
message("results at transcript level")
annoResTran$tWald$res_Wald[annoResTran$tWald$res_Wald$isoform == "PB.5435.41",]
message("results at RNA isoform level")
annoResTran$pWald$res_Wald[annoResTran$pWald$res_Wald$isoform == "PB.5435.41",]


# Ank1
redundant_iso = proteinInput$t.class.files[proteinInput$t.class.files$base_acc == "PB.9787.3", "isoform"]
IsoList <- data.frame(
  Isoform = unlist(IsoList <- list(
    subset = redundant_iso,
    protein = as.character(unique(proteinInput$peptide_orf[proteinInput$peptide_orf$isoform %in% redundant_iso,"transcript_id"]))
  )),
  Category = rep(names(IsoList), lengths(IsoList))
)
IsoList$colour <- c(rep(NA,length(IsoList$Category[IsoList$Category != "DTE"])))
pAnk1Tracks <- ggTranPlots(proteinInput$merged_peptide_gtf, proteinInput$t.class.files,isoList = c(as.character(IsoList$Isoform)), selfDf = IsoList, gene = "ANK1")
pAnk1Exp <- plot_trans_exp_individual("PB.9787.3",annoResTran$pWald$norm_counts %>% mutate(value = normalised_counts))
plot_grid(pAnk1Exp, pAnk1Tracks)

proteinInput$t.class.files[proteinInput$t.class.files$base_acc == "PB.9787.3","isoform"]
message("results at transcript level")
annoResTran$tWald$res_Wald[annoResTran$tWald$res_Wald$isoform == "PB.9787.3",]
message("results at RNA isoform level")
annoResTran$pWald$res_Wald[annoResTran$pWald$res_Wald$isoform == "PB.9787.3",]

proteinInput$t.class.files[proteinInput$t.class.files$base_acc == "PB.9787.3",] %>% 
  dplyr::select(isoform,contains("FL.")) %>% reshape2::melt(variable.name="sample",value.name="counts") %>% 
  ggplot(., aes(x = reorder(isoform, -counts), y = counts)) + geom_boxplot() + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Transcript", y = "FL reads") 



redundant_iso = proteinInput$t.class.files[proteinInput$t.class.files$associated_gene == "APOE", "isoform"]



