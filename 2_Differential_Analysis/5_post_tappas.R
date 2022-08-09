## ---------- Script -----------------
##
## Script name: 
##
## Purpose of script: Functions script for ADBDR dataset 
##
## Author: Szi Kay Leung
##
## Email: S.K.Leung@exeter.ac.uk
##
## ---------- Notes -----------------
##
## 
##   
##
## 

## ---------- Source function and config files -----------------

OUTPUT_DIR = "/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/AD_BDR/01_figures_tables"

source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR/2_Differential_Analysis/0_source_differential_functions.R")
source("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/scripts/AD_BDR/2_Differential_Analysis/bdr_differential.config.R")


## ---------- Load tappAS files -----------------
loaded <- list(
  iso = input_tappasfiles(TAPPAS_INPUT_DIR$iso),
  rna = input_tappasfiles(TAPPAS_INPUT_DIR$rna),
  iso_nf = input_tappasfiles(TAPPAS_INPUT_DIR$iso_nf),
  rna_public = input_tappasfiles(TAPPAS_INPUT_DIR$rna_public),
  rna_public_nf = input_tappasfiles(TAPPAS_INPUT_DIR$rna_public_nf)
)

# number of transcripts removed using tappAS filtering
filtered_p <- num_tappas_filter(loaded$iso$input_normalized_matrix, targeted.class.files)


## ---------- Annotate tappAS files -----------------
annotated <- list(
  iso = annotate_tappasfiles(class.files,loaded$iso$input_normalized_matrix,phenotype$iso),
  rna = annotate_tappasfiles(class.files,loaded$rna$input_normalized_matrix,phenotype$rna),
  iso_nf = annotate_tappasfiles(class.files,loaded$iso_nf$input_normalized_matrix,phenotype$iso),
  rna_public = annotate_tappasfiles(class.files,loaded$rna_public$input_normalized_matrix,phenotype$rna_public),
  rna_public_nf = annotate_tappasfiles(class.files,loaded$rna_public_nf$input_normalized_matrix,phenotype$rna_public)
)

annotated$rna_public$Norm_transcounts <- annotated$rna_public$Norm_transcounts %>% mutate(group = ifelse(group == "Control","CONTROL","CASE"))


## ---------- Gene and Transcript expression plots  -----------------
gene_exp_p <- list(
  iso = generate_plots(TargetGene,annotated$iso,"Gene","Iso-Seq Gene Expression"),
  rna = generate_plots(TargetGene,annotated$rna,"Gene","RNA-Seq Gene Expression"),
  iso_nf = generate_plots(TargetGene,annotated$iso_nf,"Gene","Iso-Seq Gene Expression"),
  rna_public = generate_plots(TargetGene,annotated$rna_public,"Gene","RNA-Seq Gene Expression"),
  rna_public_nf = generate_plots(TargetGene,annotated$rna_public_nf,"Gene","RNA-Seq Gene Expression")
)

# Expression of all transcripts 
trans_exp_p <- list(
  iso = generate_plots(TargetGene,annotated$iso,"Transcript","Iso-Seq Transcript Expression"),
  rna = generate_plots(TargetGene,annotated$rna,"Transcript","RNA-Seq Transcript Expression"),
  iso_nf = generate_plots(TargetGene,annotated$iso_nf,"Transcript","Iso-Seq Transcript Expression"),
  rna_public = generate_plots(TargetGene,annotated$rna_public,"Transcript","RNA-Seq Transcript Expression"),
  rna_public_nf = generate_plots(TargetGene,annotated$rna_public_nf,"Transcript","RNA-Seq Transcript Expression")
)


## ---------- Differential transcript expression  -----------------

difftrans <- list(
  isoseq = diff_trans_stats(loaded$iso$result_trans,targeted.class.files),
  rnaseq = diff_trans_stats(loaded$rna$result_trans,targeted.class.files),
  rnaseq_public_nf = diff_trans_stats(loaded$rna_public_nf$result_trans, targeted.class.files)
)

top20_diff_trans <- list(
  iso = generate_plots(difftrans$isoseq$transcript[1:20],annotated$iso,"Per_Transcript"),
  rna = generate_plots(difftrans$rnaseq$transcript[1:20],annotated$rna,"Per_Transcript")
)
intersect(difftrans$isoseq$transcript,difftrans$rnaseq$transcript)


## ---------- Output Plots ----------------
pdf(paste0(output_plot_dir,"/DifferentialAnalysis_Updated2.pdf"), width = 10, height = 15)
# gene expression
plot_grid(gene_exp_p$iso$APP,gene_exp_p$rna$APP,
          gene_exp_p$iso$MAPT,gene_exp_p$rna$MAPT,
          gene_exp_p$iso$FUS,gene_exp_p$rna$FUS,
          gene_exp_p$iso$SNCA,gene_exp_p$rna$SNCA,
          gene_exp_p$iso$TARDBP,gene_exp_p$rna$TARDBP, ncol = 2, 
          labels = c("A","B","C","D","E","F","G","H","I","J"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(gene_exp_p$iso$CLU,gene_exp_p$rna$CLU,
          gene_exp_p$iso$TREM2,gene_exp_p$rna$TREM2,
          gene_exp_p$iso$CD33,gene_exp_p$rna$CD33,
          gene_exp_p$iso$PTK2B,gene_exp_p$rna$PTK2B,
          gene_exp_p$iso$BIN1,gene_exp_p$rna$BIN1, ncol = 2,
          labels = c("A","B","C","D","E","F","G","H","I","J"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(gene_exp_p$iso$APOE,gene_exp_p$rna$APOE,
          gene_exp_p$iso$ABCA7,gene_exp_p$rna$ABCA7,
          gene_exp_p$iso$ABCA1,gene_exp_p$rna$ABCA1,
          gene_exp_p$iso$PICALM,gene_exp_p$rna$PICALM,
          gene_exp_p$iso$SORL1,gene_exp_p$rna$SORL1,ncol = 2,
          labels = c("A","B","C","D","E","F","G","H","I","J"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
plot_grid(gene_exp_p$iso$ANK1,gene_exp_p$rna$ANK1,
          gene_exp_p$iso$RHBDF2,gene_exp_p$rna$RHBDF2,NULL,NULL,NULL,NULL,ncol = 2,
          labels = c("A","B","C","D"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9)
# transcript expression 
plot_grid(trans_exp_p$iso$SORL1,trans_exp_p$rna$SORL1, nrow = 2)
plot_grid(trans_exp_p$iso$MAPT,trans_exp_p$rna$MAPT, nrow = 2)
plot_grid(trans_exp_p$iso$BIN1,trans_exp_p$rna$BIN1, nrow = 2)
plot_grid(trans_exp_p$iso$TARDBP,trans_exp_p$rna$TARDBP, nrow = 2)
plot_grid(trans_exp_p$iso$APP,trans_exp_p$rna$APP,nrow = 2)
plot_grid(trans_exp_p$iso$TARDBP,trans_exp_p$rna$TARDBP, nrow = 2)
plot_grid(trans_exp_p$iso$ABCA7,trans_exp_p$rna$ABCA7, nrow = 2)
plot_grid(trans_exp_p$iso$PTK2B,trans_exp_p$rna$PTK2B, nrow = 2)
plot_grid(trans_exp_p$iso$ANK1,trans_exp_p$rna$ANK1, nrow = 2)
plot_grid(trans_exp_p$iso$FYN,trans_exp_p$rna$FYN, nrow = 2)
plot_grid(trans_exp_p$iso$CLU,trans_exp_p$rna$CLU, nrow = 2)
plot_grid(trans_exp_p$iso$CD33,trans_exp_p$rna$CD33, nrow = 2)
plot_grid(trans_exp_p$iso$FUS,trans_exp_p$rna$FUS,nrow = 2)
plot_grid(trans_exp_p$iso$PICALM,trans_exp_p$rna$PICALM, nrow = 2)
plot_grid(trans_exp_p$iso$SNCA,trans_exp_p$rna$SNCA, nrow = 2)
plot_grid(trans_exp_p$iso$APOE,trans_exp_p$rna$APOE, nrow = 2)
plot_grid(trans_exp_p$iso$RHBDF2,trans_exp_p$rna$RHBDF2, nrow = 2)
plot_grid(trans_exp_p$iso$TREM2,trans_exp_p$rna$TREM2, nrow = 2)
# differential transcript expression
plot_grid(plotlist =top20_diff_trans$iso[1:10], ncol = 2, nrow = 5, 
          labels = c("A","B","C","D","E","F","G","H","I","J"),label_size = 20, label_fontfamily = "CM Roman")
plot_grid(plotlist =top20_diff_trans$iso[11:20], ncol = 2, nrow = 5,
          labels = c("K","L","M","N","O","P","Q","R","S","T"),label_size = 20, label_fontfamily = "CM Roman")
plot_grid(plotlist = top20_diff_trans$rna[1:10], ncol = 2, nrow = 5,
          labels = c("A","B","C","D","E","F","G","H","I","J"),label_size = 20, label_fontfamily = "CM Roman")
plot_grid(plotlist = top20_diff_trans$rna[11:20], ncol = 2, nrow = 5,
          labels = c("K","L","M","N","O","P","Q","R","S","T"),label_size = 20, label_fontfamily = "CM Roman")
plot_grid(filtered_p[[1]],filtered_p[[2]],NULL,NULL, labels = c("A","B"),label_size = 20, label_fontfamily = "CM Roman", scale= 0.9, rel_widths = c(0.4,0.6))
dev.off()

pdf(paste0(OUTPUT_DIR,"/DifferentialAnalysis_LargeRNASeq.pdf"), width = 10, height = 15)
for(gene in TargetGene){
  print(gene)
  print(plot_grid(trans_exp_p$iso_nf[[gene]], trans_exp_p$rna_public_nf[[gene]], 
                  ncol = 1,labels = c("A","B"), label_size = 30, label_fontfamily = "CM Roman", scale = 0.9))
}
dev.off()


## ---------- t-tests ----------------
# signficantly expressed genes
gene_stats(loaded$iso$gene_matrix, phenotype$iso)
gene_stats(loaded$rna$gene_matrix, phenotype$rna_all)

# significantly expressed transcripts 
trans_sig_ttest <- list(
  iso = run_trans_stats(annotated$iso$Norm_transcounts),
  rna = run_trans_stats(annotated$rna$Norm_transcounts)
)

pdf(paste0(output_plot_dir,"/DifferentialAnalysis_WIP.pdf"), width = 10, height = 15)
marrangeGrob(ttest_iso_plots$plots, nrow=5, ncol=2)
marrangeGrob(ttest_rna_plots$plots, nrow=5, ncol=2)
dev.off()


## ---------- RNA-Seq large public datasets (Mayo/ Rosmap) ---------------- 

plot_transexp_overtime_filter("SORL1",annotated$rna_public_nf$Norm_transcounts, loaded$rna_public_nf$input_normalized_matrix,phenotype$rna_public,"RNA-Seq")
plot_transexp_overtime_filter("SORL1",annotated$iso_nf$Norm_transcounts, loaded$iso_nf$input_normalized_matrix,phenotype$iso,"Iso-Seq")

difftrans <- c("PB.2142.35","PB.6045.790","PB.6703.1128","PB.6703.316","PB.9744.223","PB.9744.388")
difftrans_plots <- lapply(difftrans, function(x) 
  plot_grid(plot_trans_exp_individual(x, annotated$rna_public_nf$Norm_transcounts),
            plot_trans_exp_individual(x, annotated$iso_nf$Norm_transcounts))
)
difftrans_plots

tardbp_anomalies <- c("PB.94.2","PB.94.183","PB.94.239")
ank1_anomalies <- c("PB.9787.99","PB.9787.125","PB.9787.4","PB.9787.180")
for(i in ank1_anomalies){
  cat("##############",i,"#################\n")
  cat("### RNA-Seq \n")
  identify_ranked_expression(loaded$rna_public_nf$input_normalized_matrix,
                             annotated$rna_public_nf$Norm_transcounts,
                             "ANK1",phenotype$rna_public,i)
  cat("### Iso-Seq \n")
  identify_ranked_expression(loaded$iso_nf$input_normalized_matrix,
                             annotated$iso$Norm_transcounts,
                             "ANK1",phenotype$iso,i)
}

corr_filter_plot("ANK1")
corr_filter_plot("MAPT")