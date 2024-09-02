library("dplyr")
library("ggplot2")
library("ggtranscript")
library("cowplot")

collGroup <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu/AllBDRTargeted.collapsed.group.txt")
collStats <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/7_tofu/AllBDRTargeted.collapsed.read_stat.txt", header = T)

classFile <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/9_sqanti3/basic/AllBDRTargeted.collapsed_classification.filtered_lite_classification.txt", header = T)

refAllelTrancript <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals/chr17_7647112_refAlleleTranscripts.txt")
PBTrancript <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals/chr17_7647112_refAlleleTranscripts_PB.txt")
classFile[classFile$isoform %in% PBTrancript$V1,]
CCSTranscript <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals/chr17_7647112_refAlleleTranscripts_ccs.txt")[["V1"]]
CCS_c_Transcript <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals/chr17_7647112_refAlleleTranscripts_ccs_C.txt")[["V1"]]
CCS_g_Transcript <- read.table("/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/phasing/lorals/chr17_7647112_refAlleleTranscripts_ccs_G.txt")[["V1"]]

setdiff(c(as.character(CCS_c_Transcript), as.character(CCS_G_Transcript)), as.character(CCSTranscript))
setdiff(as.character(CCSTranscript), c(as.character(CCS_c_Transcript), as.character(CCS_G_Transcript)))

tcolStats <- collStats[collStats$pbid %in% PBTrancript$V1,]
tcolStats$id <- stringr::str_remove(tcolStats$id,"/ccs")

identify_reads <- function(isoform){
  
  if(isoform %in% CCS_c_Transcript){return("C")
  }else if(isoform %in% CCS_g_Transcript){return("G")
  }else{return("NA")}
  
}
tcolStats$allele <- apply(tcolStats, 1, function(x) identify_reads(x[["id"]]))
tcolStatsAllele <- tcolStats %>% group_by(pbid, allele) %>% tally() 
tcolStatsAlleleSum <- tcolStats %>% group_by(pbid) %>% tally() 
merge(tcolStatsAllele, tcolStatsAlleleSum, by = "pbid") %>%
  mutate(perc = n.x/n.y * 100) %>% 
  ggplot(., aes(x = pbid, y = perc, fill = allele)) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# obtain gtf
dirnames <- list(
  targ_iso_root = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/A_IsoSeq/",
  targ_ont_root = "/lustre/projects/Research_Project-MRC148213/lsl693/AD_BDR/D_ONT/")
gtf <- list(
  iso = rtracklayer::import(paste0(dirnames$targ_iso_root, "9_sqanti3/full/AllBDRTargeted.collapsed_classification.filtered_lite.gtf")),
  ont_merged = rtracklayer::import(paste0(dirnames$targ_ont_root, "5_cupcake/7_sqanti3/ontBDR_collapsed.filtered_counts_filtered.gtf"))
)
gtf <- lapply(gtf, function(x) as.data.frame(x))

plot_allele_perTranscript <- function(inputgtf, transcript, refValue, altValue){
  dat <- inputgtf[inputgtf$transcript_id == transcript,]
  
  dat <- dat[dat$type == "exon",] %>% mutate(value = refValue)
  
  inter_exons <- dat %>%
    mutate(
      next_start = lead(start) - 1
    ) %>%
    filter(!is.na(next_start)) %>%
    transmute(
      seqnames = seqnames,
      start = end + 1,
      end = next_start,
      width = end - start + 1,
      strand = strand,
      source = source,
      type = "inter_exon",
      score = score,
      phase = phase,
      gene_id = gene_id,
      transcript_id = transcript_id
    )
  
  # Combine exon and inter-exon data frames
  ref_data <- bind_rows(dat, inter_exons) %>% arrange(start) %>% mutate(transcript = transcript)
  
  ref_lines <- ref_data %>%
    select(seqnames, start, end) %>%
    tidyr::pivot_longer(cols = c(start, end), names_to = "position_type", values_to = "position") %>%
    mutate(value = ifelse(position_type == "start", refValue, 0.01),
           transcript = transcript)
  
  ref_data <- ref_data %>% mutate(value = ifelse(is.na(value),0.01,refValue), transcript = transcript)
  alt_data <- ref_data %>% mutate(value = ifelse(value == 0.01,0.01,altValue), transcript = transcript)
  alt_lines <- ref_lines %>% mutate(value = ifelse(value == 0.01,0.01,altValue))
  
  output <- list(ref_data, ref_lines, alt_data, alt_lines)
  names(output) <- c("ref","ref_lines","alt", "alt_lines")
  return(output)
}

output <- plot_allele_perTranscript(gtf$iso, "PB.4702.13", 7, 1)
# Plot the data using ggplot2
ggplot() +
  geom_segment(data = output$ref, aes(x = start, xend = end, y = value, yend = value), linetype = "solid", colour = "blue") +
  geom_segment(data = output$ref_lines, aes(x = position, xend = position, y = 0, yend = value), linetype = "dotted", colour = "blue") +
  geom_segment(data = output$alt_lines, aes(x = position, xend = position, y = 0, yend = value), linetype = "dotted", colour = "orange") +
  geom_segment(data = output$alt, aes(x = start, xend = end, y = value, yend = value), linetype = "solid", colour = "orange") +
  labs(x = "Position", y = "Value") +
  theme_classic() +
  theme(axis.text.y = element_text(angle=90)) 

p2 <- dat %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range() + labs(x = NULL, y = "") +
  geom_intron(data = to_intron(dat, "transcript_id"),
              aes(strand = strand)) +
  theme_classic() +
  theme(axis.text.y = element_text(angle=90)) 

plot_grid(p1,p2,nrow=2)

tcolStatsAllele[tcolStatsAllele$pbid == "PB.4702.13",]
# remove NA alleles
tcolStatsAllele <- tcolStatsAllele[tcolStatsAllele$allele != "NA",]
tcolStatsAlleleW <- tidyr::spread(tcolStatsAllele, "allele", "n") 
tcolStatsAlleleW[is.na(tcolStatsAlleleW)] <- 0

tcolStatsAlleleW <- as.data.frame(tcolStatsAlleleW) %>% mutate(totalReads = tcolStatsAlleleW$C + tcolStatsAlleleW$G)
tcolStatsAlleleW <- tcolStatsAlleleW[tcolStatsAlleleW$totalReads > 20,]
plotOutput <- apply(tcolStatsAlleleW, 1, function(x) plot_allele_perTranscript(gtf$iso, x[["pbid"]], x[["C"]],x[["G"]]))
names(plotOutput) <- tcolStatsAlleleW$pbid

refs <- lapply(plotOutput, function(x) x$ref)
combined_refs <- do.call(rbind, refs)
combined_refs <- combined_refs %>% mutate(transcript = factor(transcript, levels = rev(tcolStatsAlleleW$pbid)))

ref_lines <- lapply(plotOutput, function(x) x$ref_lines)
combined_refs_lines <- do.call(rbind, ref_lines)

alts <- lapply(plotOutput, function(x) x$alt)
combined_alts <- do.call(rbind, alts)

alt_lines <- lapply(plotOutput, function(x) x$alt_lines)
combined_alt_lines <- do.call(rbind, alt_lines)

p1 <- ggplot() +
  geom_segment(data = combined_refs, aes(x = start, xend = end, y = as.numeric(value), yend = as.numeric(value)), linetype = "solid", colour = "blue") +
  geom_segment(data = combined_refs_lines, aes(x = position, xend = position, y = 0, yend = as.numeric(value)), linetype = "dotted", colour = "blue") +
  geom_segment(data = combined_alt_lines, aes(x = position, xend = position, y = 0, yend = as.numeric(value)), linetype = "dotted", colour = "orange") +
  geom_segment(data =  combined_alts, aes(x = start, xend = end, y = as.numeric(value), yend = as.numeric(value)), linetype = "solid", colour = "orange") +
  facet_grid(transcript ~ .) +
  labs(x = "Position", y = "Value") +
  theme_classic() +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 44))# +
  #theme(axis.text.y = element_text(angle=90)) 

datGtf <- gtf$iso %>%
  filter(transcript_id %in% tcolStatsAlleleW$pbid & type == "exon") %>%
  mutate(transcript_id = factor(transcript_id, levels = rev(tcolStatsAlleleW$pbid)))
p2 <- datGtf %>%
  ggplot(aes(
    xstart = start,
    xend = end,
    y = transcript_id
  )) +
  geom_range() + labs(x = NULL, y = "") +
  geom_intron(data = to_intron(datGtf, "transcript_id"),
              aes(strand = strand)) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text on the left
    axis.ticks.y = element_blank(),  # Remove y-axis ticks on the left
    axis.line.y = element_blank(),  # Remove y-axis line on the left
    axis.text.y.left = element_text(color = "black")  # Add y-axis text on the right
  ) #+
  #scale_y_discrete(
  #  position = "right",  # Move y-axis to the right
  #)

plot_grid(p1,p2, nrow=2)
