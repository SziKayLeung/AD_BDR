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
    mutate(value = ifelse(position_type == "start", refValue, NA),
           transcript = transcript)
  
  # Append the final end position to ref_lines
  last_end <- data.frame(position = max(dat$end), value = refValue, transcript = transcript)
  ref_lines <- bind_rows(ref_lines, last_end)
  
  ref_data <- ref_data %>% mutate(value = ifelse(is.na(value), NA, refValue), transcript = transcript)
  alt_data <- ref_data %>% mutate(value = ifelse(is.na(value), NA, altValue), transcript = transcript)
  alt_lines <- ref_lines %>% mutate(value = ifelse(is.na(value), NA, altValue))
  
  # Append the final end position to alt_lines
  last_end_alt <- data.frame(position = max(dat$end), value = altValue, transcript = transcript)
  alt_lines <- bind_rows(alt_lines, last_end_alt)
  
  output <- list(ref_data, ref_lines, alt_data, alt_lines)
  names(output) <- c("ref", "ref_lines", "alt", "alt_lines")
  return(output)
}