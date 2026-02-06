#' Collapse events per sample and chromosome
#'
#' This function collapses genomic events per individual, chromosome, and group,
#' summarising the number of events, total Mendelian errors, the total span size,
#' and a string listing all merged event coordinates.
#'
#' @details
#' When values are available in the ratio columns (`ratio_proband`, `ratio_mother`, `ratio_father`), 
#' weighted mean ratios are computed across the collapsed events.
#' The weighted mean ratio is calculated as:
#'
#' \deqn{\frac{\sum_{i} r_i \times N_i}{\sum_{i} N_i}}
#'
#' where \eqn{r_i} is the ratio of each individual event and \eqn{N_i} the number
#' of SNPs in that event.
#' If all values in a ratio column are `NA`, the corresponding collapsed value will be `NA_real_`.
#' 
#' @param subset_df A data.frame containing event-level data with columns:
#' \itemize{
#'   \item ID – Sample identifier.
#'   \item chromosome – Chromosome name.
#'   \item start – Start position of the event.
#'   \item end – End position of the event.
#'   \item group – Event group/class.
#'   \item n_snps – Number of SNPs in the event.
#'   \item n_mendelian_error – Number of Mendelian errors in the event.
#'   \item ratio_proband, ratio_mother, ratio_father – depth-ratio metrics.
#' }
#'   
#' @param min_ME Minimum number of Mendelian errors required to retain an event
#'   before collapsing (default: 2).
#'   
#' @param min_size Minimum genomic span size required to retain an event
#'   before collapsing, in base pairs (default: 500e3).
#'   
#' @return A data.frame with collapsed events and columns:
#' \itemize{
#'   \item ID – Sample identifier.
#'   \item chromosome – Chromosome name.
#'   \item start, end – Genomic span of the collapsed block.
#'   \item group – HMM state of the block.
#'   \item n_events – Number of events collapsed.
#'   \item total_mendelian_error – Sum of Mendelian errors across events.
#'   \item total_size – Total genomic span size covered by events.
#'   \item total_snps – Total SNPs in the overlapping events.
#'   \item prop_covered – Proportion of the region covered by events.
#'   \item ratio_proband, ratio_mother, ratio_father – Weighted mean ratios across the collapsed events.
#'         If all values are NA, the column will contain NA_real_.
#'   \item collapsed_events – Comma-separated list of collapsed events.
#' }
#'
#' @export
#' @examples
#' all_events <- data.frame(
#' ID = c("S1", "S1", "S1", "S2", "S2"),
#' chromosome = c("1", "1", "1", "2", "2"),
#' start = c(100e4, 200e4, 300e4, 500e4, 600e4),
#' end   = c(160e4, 260e4, 360e4, 560e4, 700e4),
#' group = c("iso_mat", "iso_mat", "het_pat", "iso_mat", "iso_mat"),
#' n_mendelian_error = c(5, 10, 2, 50, 30),
#' stringsAsFactors = FALSE
#' )
#' out <- collapseEvents(all_events)
#' 
collapseEvents <- function(subset_df, min_ME = 2, min_size = 500e3) {
  # Create event string and size per row
  subset_df$event_string <- paste0(
    subset_df$chromosome, ":", subset_df$start, "-", subset_df$end
  )
  subset_df$event_size <- subset_df$end - subset_df$start
  
  # Filter events based on quality thresholds:
  # keep only those with ≥ min_ME Mendelian errors and ≥ min_size genomic span
  subset_df <- subset_df[
    subset_df$n_mendelian_error > min_ME &
      subset_df$event_size > min_size,
  ]
  
  # Return empty output if no events remain after filtering
  if (nrow(subset_df) == 0) {
    return(data.frame(
      ID = character(),
      chromosome = character(),
      start = numeric(),
      end = numeric(),
      group = character(),
      n_events = numeric(),
      total_mendelian_error = numeric(),
      total_size = numeric(),
      total_snps = numeric(),
      prop_covered = numeric(),
      ratio_father  = numeric(),
      ratio_mother  = numeric(),
      ratio_proband = numeric(),
      collapsed_events = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Create grouping key
  subset_df$group_key <- paste(subset_df$ID, subset_df$chromosome, subset_df$group, sep = "_")
  
  # Split by key
  splitted <- split(subset_df, subset_df$group_key)
  
  # Collapse manually
  collapsed_list <- lapply(splitted, function(df) {
    
    region_min <- min(df$start, na.rm = TRUE)
    region_max <- max(df$end,   na.rm = TRUE)
    region_span <- region_max - region_min
    
    event_sizes <- df$event_size
    snps        <- df$n_snps
 
    w <- snps
    
    ratio_father_val  <- if(all(is.na(df$ratio_father)))  NA_real_ else stats::weighted.mean(df$ratio_father,  w, na.rm = TRUE)
    ratio_mother_val  <- if(all(is.na(df$ratio_mother)))  NA_real_ else stats::weighted.mean(df$ratio_mother,  w, na.rm = TRUE)
    ratio_proband_val <- if(all(is.na(df$ratio_proband))) NA_real_ else stats::weighted.mean(df$ratio_proband, w, na.rm = TRUE)
    
    out <- data.frame(
      ID = df$ID[1],
      chromosome = df$chromosome[1],
      start = region_min,
      end = region_max,
      group = df$group[1],
      n_events = nrow(df),
      total_mendelian_error = sum(df$n_mendelian_error, na.rm = TRUE),
      total_size = region_span,
      total_snps = sum(snps, na.rm = TRUE),
      prop_covered = sum(event_sizes, na.rm = TRUE) / region_span,
      ratio_father  = ratio_father_val,
      ratio_mother  = ratio_mother_val,
      ratio_proband = ratio_proband_val,
      collapsed_events = paste(df$event_string, collapse = ","),
      stringsAsFactors = FALSE
    )
    
    return(out)
    
  })
  
  collapsed_events <- do.call(rbind, collapsed_list)
  rownames(collapsed_events) <- NULL
  
  return(collapsed_events)
}
