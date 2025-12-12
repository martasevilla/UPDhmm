#' Identify recurrent genomic regions across samples
#'
#' This function finds recurrent genomic regions across samples based on overlapping intervals.
#' Regions exceeding the Mendelian error threshold are excluded.
#'
#' @param df A data.frame with columns:
#'   - ID: Sample identifier.
#'   - chromosome: Chromosome name.
#'   - start: Start coordinate of the region.
#'   - end: End coordinate of the region.
#'   - n_mendelian_error: Number of Mendelian errors in the region.
#' @param error_threshold Numeric, default = 100.
#'   Maximum number of Mendelian errors allowed for a region to be considered.
#' @param min_support Integer, default = 3.
#'   Minimum number of unique samples required to call a region recurrent.
#' @param ID_col Character string indicating the column name containing
#'   sample identifiers. Default is `"ID"`.
#'
#' @return A GRanges object containing the recurrent regions that meet
#'   the minimum support threshold.
#' @export
#' @examples
#' df <- data.frame(
#' ID = c("S1", "S2", "S3", "S3"),
#' chromosome = c("chr1", "chr1", "chr1", "chr1"),
#' start = c(100, 120, 500, 510),
#' end = c(150, 170, 550, 560),
#' n_mendelian_error = c(10, 20, 5, 5)
#' )
#' identifyRecurrentRegions(df, ID_col = "ID", error_threshold = 50, min_support = 2)

identifyRecurrentRegions <- function(df,
                                     ID_col = "ID",
                                     error_threshold = 100,
                                     min_support = 3) {
  #---------------------------------------------------------------
  # 1. Validate inputs
  #---------------------------------------------------------------
  if (!ID_col %in% names(df)) {
    stop(sprintf("Column '%s' not found in input dataframe.", ID_col))
  }
  
  if (all(c("start", "end") %in% colnames(df))) {
    start_col <- "start"
    end_col   <- "end"
  } else if (all(c("min_start", "max_end") %in% colnames(df))) {
    start_col <- "min_start"
    end_col   <- "max_end"
  } else {
    stop("Input must contain either (start, end) or (min_start, max_end) columns.")
  }
  
  err_col <- if ("n_mendelian_error" %in% names(df)) "n_mendelian_error" else
    if ("total_mendelian_error" %in% names(df)) "total_mendelian_error" else
      stop("Input must contain either n_mendelian_error or total_mendelian_error column.")
  
  #---------------------------------------------------------------
  # 2. Convert to GRanges
  #---------------------------------------------------------------
  gr <- GenomicRanges::GRanges(
    seqnames = df$chromosome,
    ranges = IRanges::IRanges(
      start = as.numeric(df[[start_col]]),
      end = as.numeric(df[[end_col]])
    ),
    ID = df[[ID_col]],
    n_mendelian_error = df[[err_col]]
  )
  
  #---------------------------------------------------------------
  # 3. Filter and merge overlapping regions
  #---------------------------------------------------------------
  gr_filtered <- gr[S4Vectors::mcols(gr)$n_mendelian_error < error_threshold]
  
  if (length(gr_filtered) == 0)
    return(GenomicRanges::GRanges())
  
  reduced_gr <- GenomicRanges::reduce(gr_filtered)
  
  #---------------------------------------------------------------
  # 4. Count distinct samples per merged region
  #---------------------------------------------------------------
  hits <- GenomicRanges::findOverlaps(reduced_gr, gr_filtered)
  hits_df <- data.frame(
    region_index = S4Vectors::queryHits(hits),
    ID = S4Vectors::mcols(gr_filtered)$ID[S4Vectors::subjectHits(hits)],
    stringsAsFactors = FALSE
  )
  
  region_support <- stats::aggregate(ID ~ region_index, data = hits_df,
                                     FUN = function(x) length(unique(x)))
  
  #---------------------------------------------------------------
  # 5. Keep only recurrent regions
  #---------------------------------------------------------------
  recurrent_regions <- reduced_gr[region_support$ID >= min_support]
  recurrent_regions$n_samples <- region_support$ID[region_support$ID >= min_support]
  
  return(recurrent_regions)
}
