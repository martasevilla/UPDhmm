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
#' @param largeCollapsedVcf Optional CollapsedVCF object. If provided, depth ratios
#'   per sample (proband, mother, father) will be calculated for the collapsed events.
#'
#' @param field_DP Optional character string specifying which VCF FORMAT field to use 
#'   for depth metrics (e.g., "DP" or "AD"). Default = NULL, which will try "DP" first 
#'   and then "AD".
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
#'   collapsed_events: Comma-separated list of collapsed events
#'   min_start, max_end: Genomic span of collapsed block
#'   ratio_proband, ratio_mother, ratio_father: Normalized depth ratios per sample
#'     (only present if `largeCollapsedVcf` is provided).
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
collapseEvents <- function(subset_df, largeCollapsedVcf = NULL, field_DP = NULL, min_ME = 2, min_size = 500e3) {
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
    return(data.frame(
      ID = character(),
      seqnames = character(),
      group = character(),
      n_events = numeric(),
      total_mendelian_error = numeric(),
      total_size = numeric(),
      collapsed_events = character(),
      min_start = numeric(),
      max_end = numeric(),
      ratio_proband = numeric(),
      ratio_mother = numeric(),
      ratio_father = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
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
  
  # Compute depth ratios if VCF provided
  if (!is.null(largeCollapsedVcf)) {
    
    geno_all <- VariantAnnotation::geno(largeCollapsedVcf)
    
    dp_field <- if (!is.null(field_DP) && field_DP %in% names(geno_all)) { field_DP }
    else if ("DP" %in% names(geno_all)) { "DP" }
    else if ("AD" %in% names(geno_all)) { "AD" }
    else stop("No DP/AD field found in VCF")
    
    if (dp_field == "AD") {
      ad <- geno_all$AD
      depth_matrix <- sapply(seq_len(ncol(ad)), function(j) {
        vapply(ad[, j], function(x) if (all(is.na(x))) NA else sum(x, na.rm = TRUE), numeric(1))
      })
      depth_matrix <- as.matrix(depth_matrix)
      colnames(depth_matrix) <- colnames(ad)
    } else {
      depth_matrix <- as.matrix(geno_all[[dp_field]])
    }
    
    total_mean <- computeTrioTotals(largeCollapsedVcf, field_DP = field_DP)

    gr_events <- GenomicRanges::GRanges(
      seqnames = collapsed_events$seqnames,
      ranges = IRanges(collapsed_events$min_start, collapsed_events$max_end)
    )
    
    hits <- GenomicRanges::findOverlaps(gr_events, rowRanges(largeCollapsedVcf))
    idx_list <- split(S4Vectors::subjectHits(hits), S4Vectors::queryHits(hits))
    
    compute_ratio_event <- function(vcf_idx) {
      if (length(vcf_idx) == 0) return(rep(NA, length(total_mean)))
      dp <- depth_matrix[vcf_idx, , drop = FALSE]
      block_mean <- colMeans(dp, na.rm = TRUE)
      block_mean / total_mean
    }
    
    ratio_list <- lapply(idx_list, compute_ratio_event)
    
    # Fill NA for non-overlapping events
    fill_na <- rep(list(rep(NA, 3)), nrow(collapsed_events))
    fill_na[as.integer(names(ratio_list))] <- ratio_list
    
    collapsed_events$ratio_proband <- sapply(fill_na, `[`, 1)
    collapsed_events$ratio_mother  <- sapply(fill_na, `[`, 2)
    collapsed_events$ratio_father  <- sapply(fill_na, `[`, 3)
  }
  
  chr_num <- suppressWarnings(as.numeric(collapsed_events$seqnames))
  collapsed_events <- collapsed_events[order(collapsed_events$ID, chr_num), ]
  rownames(collapsed_events) <- NULL
  
  return(collapsed_events)
}
