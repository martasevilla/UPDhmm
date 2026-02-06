#' Annotate regions as recurrent or non-recurrent
#'
#' Given a results data.frame and a set of recurrent genomic regions,
#' this function labels each row as "Yes" (recurrent) or "No" (non-recurrent)
#' based on overlaps with a set of recurrent regions.
#'
#' @param df Data.frame with region coordinates and sample IDs.
#' @param recurrent_regions GRanges object from identifyRecurrentRegions().
#' @param min_overlap Numeric between 0 and 1, default = 0.7.
#'    Minimum fraction of the input region that must overlap a recurrent region to be annotated as recurrent.
#'   
#' @details
#' A region is marked as recurrent if the fraction of its length overlapping a recurrent region 
#' is at least `min_overlap` (default 0.7).
#' 
#' @return The same data.frame with two added columns:
#'   - Recurrent: "Yes" or "No"
#'   - n_samples: Number of supporting samples (if recurrent)
#'   
#' @export
#' @examples
#' input <- data.frame(
#'     ID = c("S1", "S2", "S3", "S4"),
#'     chromosome = c("chr1", "chr1", "chr1", "chr2"),
#'     start = c(100, 150, 500, 100),
#'     end = c(150, 200, 550, 150),
#'     n_mendelian_error = c(10, 20, 5, 200)
#' )
#'
#'recurrent_gr <- GenomicRanges::GRanges(
#'  seqnames = "chr1",
#'  ranges = IRanges::IRanges(
#'    start = 100,
#'    end = 170
#'  ),
#'  n_samples = 2
#')
#' markRecurrentRegions(input, recurrent_gr)
#' 
markRecurrentRegions <- function(df, recurrent_regions, min_overlap = 0.7) {
  
  # --- Detect coordinate columns automatically ---
  if (!all(c("start", "end") %in% names(df))) {
    stop("Input must contain columns 'start' and 'end'.")
  }
  
  # --- Initialize default output columns ---
  df$Recurrent <- "No"
  df$n_samples <- 1
  
  # --- If no recurrent regions are provided, return input unchanged ---
  if (length(recurrent_regions) == 0) {
    return(df)
  }
  
  # --- Convert input table to GRanges ---
  gr <- GenomicRanges::GRanges(
    seqnames = df$chromosome,
    ranges   = IRanges::IRanges(df$start, df$end),
    ID       = df$ID
  )
  
  # --- Identify overlaps with recurrent regions ---
  overlaps <- GenomicRanges::findOverlaps(gr, recurrent_regions)
  
  # --- If no overlaps are found, return input unchanged ---
  if (length(overlaps) == 0) {
    return(df)
  }
  
  qh <- S4Vectors::queryHits(overlaps)
  sh <- S4Vectors::subjectHits(overlaps)
  
  # --- Compute overlap fraction for each overlapping event ---
  
  # Compute intersection between input regions and recurrent regions
  ov <- GenomicRanges::pintersect(gr[qh], recurrent_regions[sh])
  
  # Fraction of the input event covered by the recurrent region
  frac_covered <- IRanges::width(ov) / IRanges::width(gr[qh])
  
  # Retain only events exceeding the minimum overlap threshold
  keep <- frac_covered >= min_overlap
  
  # --- Annotate recurrent events ---
  if (any(keep)) {
    df$Recurrent[qh[keep]] <- "Yes"
    df$n_samples[qh[keep]] <-
      recurrent_regions$n_samples[sh[keep]]
  }
  
  rownames(df) <- NULL
  return(df)
}
