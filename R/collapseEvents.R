#' Collapse events per sample and chromosome
#'
#' This function collapses genomic events per individual, chromosome, and group,
#' summarising the number of events, total Mendelian errors, the total span size,
#' and a string listing all merged event coordinates.
#'
#' @param subset_df A data.frame containing event-level data with columns:
#'   ID: Sample identifier.
#'   seqnames: Chromosome name.
#'   start: Start position of the event.
#'   end: End position of the event.
#'   group: Event group/class.
#'   n_snps: Number of SNPs in the event. 
#'   n_mendelian_error: Number of Mendelian errors in the event.
#'   ratio_proband, ratio_mother, ratio_father: (optional) depth-ratio metrics
#'   
#' @param min_ME Minimum number of Mendelian errors required to retain an event
#'   before collapsing (default: 2).
#'   
#' @param min_size Minimum genomic span size required to retain an event
#'   before collapsing, in base pairs (default: 500e3).
#'   
#' @return A data.frame with collapsed events and columns:
#'  ID, seqnames, group
#'  n_events: Number of events collapsed
#'  total_mendelian_error: Sum of Mendelian errors across events
#'  total_size: Total genomic span size covered by events
#'  collapsed_events: Comma-separated list of collapsed events
#'  min_start, max_end: Genomic span of collapsed block
#'  total_snps: Total SNPs in the overlapping events
#'  prop_covered: Proportion of the region covered by events
#'  ratio_proband, ratio_mother, ratio_father: Weighted mean ratios across the collapsed events,
#'    calculated as \(\sum_i r_i \cdot N_i / \sum_i N_i\), where \(r_i\) is the ratio of each
#'    individual event and \(N_i\) is the number of SNPs in that event. Returned only if
#'    the input contains all three ratio columns
#'
#' @export
#' @examples
#' all_events <- data.frame(
#' ID = c("S1", "S1", "S1", "S2", "S2"),
#' seqnames = c("1", "1", "1", "2", "2"),
#' start = c(100, 150, 300, 500, 550),
#' end = c(120, 180, 320, 520, 580),
#' group = c("iso_mat", "iso_mat", "het_pat", "iso_mat", "iso_mat"),
#' n_mendelian_error = c(5, 10, 2, 50, 30),
#' stringsAsFactors = FALSE
#' )
#' out <- collapseEvents(all_events)
collapseEvents <- function(subset_df, min_ME = 2, min_size = 500e3) {
  # Create event string and size per row
  subset_df$event_string <- paste0(
    subset_df$seqnames, ":", subset_df$start, "-", subset_df$end
  )
  subset_df$event_size <- subset_df$end - subset_df$start
  
  # Filter events based on quality thresholds:
  # keep only those with ≥ min_ME Mendelian errors and ≥ min_size genomic span
  subset_df <- subset_df[
    subset_df$n_mendelian_error >= min_ME &
      subset_df$event_size >= min_size,
  ]
  
  # Return empty output if no events remain after filtering
  if (nrow(subset_df) == 0) {
    
    base_cols <- list(
      ID = character(),
      seqnames = character(),
      group = character(),
      n_events = numeric(),
      total_mendelian_error = numeric(),
      total_size = numeric(),
      collapsed_events = character(),
      min_start = numeric(),
      max_end = numeric(),
      total_snps = numeric(),
      prop_covered = numeric()
    )
    
    # Add ratio columns ONLY if present in input
    if (all(c("ratio_proband", "ratio_mother", "ratio_father") %in% colnames(subset_df))) {
      base_cols$ratio_proband <- numeric()
      base_cols$ratio_mother  <- numeric()
      base_cols$ratio_father  <- numeric()
    }
    
    return(as.data.frame(base_cols, stringsAsFactors = FALSE))
  }
  
  
  # Create grouping key
  subset_df$group_key <- paste(subset_df$ID, subset_df$seqnames, subset_df$group, sep = "_")
  
  # Split by key
  splitted <- split(subset_df, subset_df$group_key)
  
  # Collapse manually
  collapsed_list <- lapply(splitted, function(df) {
    
    region_min <- min(df$start, na.rm = TRUE)
    region_max <- max(df$end,   na.rm = TRUE)
    region_span <- region_max - region_min
    
    event_sizes <- df$event_size
    snps        <- df$n_snps
    
    out <- data.frame(
      ID = df$ID[1],
      seqnames = df$seqnames[1],
      group = df$group[1],
      n_events = nrow(df),
      total_mendelian_error = sum(df$n_mendelian_error, na.rm = TRUE),
      total_size = sum(event_sizes, na.rm = TRUE),
      collapsed_events = paste(df$event_string, collapse = ","),
      min_start = region_min,
      max_end = region_max,
      total_snps = sum(snps, na.rm = TRUE),
      prop_covered = sum(event_sizes, na.rm = TRUE) / region_span,
      stringsAsFactors = FALSE
    )
    
    if (all(c("ratio_proband", "ratio_mother", "ratio_father") %in% colnames(df))) {
      w <- snps
      out$ratio_proband <- weighted.mean(df$ratio_proband, w, na.rm = TRUE)
      out$ratio_mother  <- weighted.mean(df$ratio_mother,  w, na.rm = TRUE)
      out$ratio_father  <- weighted.mean(df$ratio_father,  w, na.rm = TRUE)
    }
    
    return(out)
    
  })
  
  collapsed_events <- do.call(rbind, collapsed_list)
  
  chr_num <- as.numeric(collapsed_events$seqnames)
  collapsed_events <- collapsed_events[order(collapsed_events$ID, chr_num), ]
  rownames(collapsed_events) <- NULL
  
  return(collapsed_events)
}
