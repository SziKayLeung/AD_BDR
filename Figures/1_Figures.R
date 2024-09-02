## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##



## ---------- Source function and config files -----------------

suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("tibble"))
suppressMessages(library("stringr"))

LOGEN_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/LOGen/"
SC_ROOT = "/lustre/projects/Research_Project-MRC148213/lsl693/scripts/AD_BDR/Figures/"
source(paste0(SC_ROOT, "BDR_config.R"))
output_dir = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/01_figures_tables"
source(paste0(LOGEN_ROOT, "/aesthetics_basics_plots/pthemes.R"))
source(paste0(LOGEN_ROOT, "differential_analysis/plot_transcript_level.R"))
sapply(list.files(path = paste0(LOGEN_ROOT,"target_gene_annotation"), pattern="*summarise*", full = T), source,.GlobalEnv)


for(i in annoResTran){
  print(paste0("Number of singificant sites:", nrow(i[["anno_res"]])))
}

annoResTran$ontMergedBraak$anno_res %>% group_by(associated_gene) %>% tally()


setorderAD = c("Control","Intermediate_1","Intermediate_2","AD")
setorderAD = c("Control","Intermediate_1","Intermediate_2","AD")
setorderBraak = c(0,seq(1:6))
plot_transexp_overtime("SNCA",annoResTran$B2Wald$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.70995.7653"), setorder=setorderAD)
plot_transexp_overtime("APP",annoResTran$B2WaldBraak$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.57985.21"),setorder=setorderBraak)
a = plot_transexp_overtime("BIN1",annoResTran$ontMergedBraakSum$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.50706.8"),setorder=setorderBraak)
b = plot_transexp_overtime("BIN1",annoResTran$B1WaldBraak$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.50706.8"),setorder=setorderBraak)
c = plot_transexp_overtime("BIN1",annoResTran$B2WaldBraak$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.50706.8"),setorder=setorderBraak)
d = plot_transexp_overtime("BIN1",annoResTran$ontMergedBraak$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.50706.8"),setorder=setorderBraak)

plot_grid(a, d, b, c, nrow = 1, labels = c("TotalSum","TotalSep","B1","B2"))

plot_transexp_overtime("APP",TargetedDESeq$ontResTranAnno$ontMergedBraakSum$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c("PB.57985.21"),setorder=setorderBraak)
Trem2_p <- list(
  dendro = plot_dendro_Tgene(dirnames$targ_anno, "TREM2"),
  tracks = NULL,
  ONTGeneExp  = generate_plots("Trem2",annotated$targ_ont,"Gene_time","")$Trem2,
  ONTTransExp = generate_plots("Trem2",annotated$targ_ont,"Transcript_Iso_Trajectory"," ",tappassigtrans$targ_ont$TargetedOnt_Transexp)$Trem2,
  IF = plot_grid(generate_plots("Trem2",loaded$targ_ont,"IF_time_series",phenotype$targ_ont,"ONT_Targeted")[[1]]),
  pheat = draw_heatmap_gene("Trem2", class.files$targ_ont, annotated$targ_ont$Norm_transcounts)$gtable
)

plot_dendro_Tgene(dirnames$targ_anno, "TREM2")
plot_dendro_Tgene(dirnames$targ_anno, "CD33")
plot_dendro_Tgene(dirnames$targ_anno, "APP")

merge(annoResTran$ontMergedBraak$norm_counts, phenotype$batch_ont, by = "sample") %>% 
  filter(isoform == "PB.50706.24") %>% 
  ggplot(., aes(x = as.factor(BraakTangle_numeric), y = normalised_counts)) + geom_boxplot(outlier.shape = NA) + geom_point() + 
  facet_grid(~APOE_geno, scales = "free") 



pdf("DifferentialAnalysis2.pdf")
for(i in 1:nrow(TargetedDESeq$ontResTranAnno$ontMergedBraak$anno_res)){
  gene = TargetedDESeq$ontResTranAnno$ontMergedBraak$anno_res$associated_gene[i]
  iso = TargetedDESeq$ontResTranAnno$ontMergedBraak$anno_res$isoform[i]
  print(plot_transexp_overtime(gene,annoResTran$ontMergedBraak$norm_counts,design="multiple_case_control",show="specific",isoSpecific=c(iso),setorder=setorderBraak))
}
dev.off()

p <- annoResTran$ontMerged$norm_counts %>% filter(isoform == iso) %>% 
  filter(group %in% c(0,6)) %>%
  ggplot(., aes(x = reorder(sample, normalised_counts), y = normalised_counts, fill = group)) +
  geom_bar(stat = "identity") +
  facet_grid(~group, scales="free") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p)
