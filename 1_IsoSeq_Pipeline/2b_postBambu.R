library("stringr")

dat <- data.table::fread("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/A_IsoSeq/7_bambu/CPM_transcript.txt")

# Specify the list of gene symbols
gene_table <- data.frame(
  GeneName = c("ABCA1", "SORL1", "MAPT", "BIN1", "TARDBP", "APP", "ABCA7", "PTK2B",
               "ANK1", "FYN", "CLU", "CD33", "FUS", "PICALM", "SNCA", "APOE", 
               "TRPA1", "RHBDF2", "TREM2", "VGF"),
  EnsemblID = c("ENSG00000165029", "ENSG00000102082", "ENSG00000186868", "ENSG00000160315", "ENSG00000120948",
                "ENSG00000142192", "ENSG00000170396", "ENSG00000149248", "ENSG00000145675", "ENSG00000102882",
                "ENSG00000120885", "ENSG00000166710", "ENSG00000090339", "ENSG00000188157", "ENSG00000145335",
                "ENSG00000130203", "ENSG00000196628", "ENSG00000147669", "ENSG00000183148", "ENSG00000140374")
)

dat$GENEID1 <- sapply(dat$GENEID, function(x) word(x,c(1),sep=fixed(".")))
targeted <- dat %>% filter(GENEID1 %in% gene_table$EnsemblID) %>% left_join(., gene_table, by = c("GENEID1" = "EnsemblID"))
