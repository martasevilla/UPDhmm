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
#'   n_mendelian_error: Number of Mendelian errors in the event.
#'
#' @return A data.frame with collapsed events and columns:
#'  ID, seqnames, group
#'  n_events: Number of events collapsed
#'  total_mendelian_error: Sum of Mendelian errors across events
#'  total_size: Total genomic span size covered by events
#'   collapsed_events: Comma-separated list of collapsed events
#'   min_start, max_end: Genomic span of collapsed block
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
collapseEvents <- function(subset_df) {
  # Create event string and size per row
  subset_df$event_string <- paste0(
    subset_df$seqnames, ":", subset_df$start, "-", subset_df$end
  )
  subset_df$event_size <- subset_df$end - subset_df$start
  
  # Create grouping key
  subset_df$group_key <- paste(subset_df$ID, subset_df$seqnames, subset_df$group, sep = "_")
  
  # Split by key
  splitted <- split(subset_df, subset_df$group_key)
  
  # Collapse manually
  collapsed_list <- lapply(splitted, function(df) {
    data.frame(
      ID = df$ID[1],
      seqnames = df$seqnames[1],
      group = df$group[1],
      n_events = nrow(df),
      total_mendelian_error = sum(df$n_mendelian_error, na.rm = TRUE),
      total_size = sum(df$event_size, na.rm = TRUE),
      collapsed_events = paste(df$event_string, collapse = ","),
      min_start = min(df$start, na.rm = TRUE),
      max_end = max(df$end, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  
  collapsed_events <- do.call(rbind, collapsed_list)
  rownames(collapsed_events) <- NULL
  return(collapsed_events)
}
