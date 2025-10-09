#' Annotate regions as recurrent or non-recurrent
#'
#' Given a results data.frame and a set of recurrent genomic regions,
#' this function labels each row as "Yes" (recurrent) or "No".
#'
#' @param df Data.frame with region coordinates and sample IDs.
#' @param recurrent_regions GRanges object from identifyRecurrentRegions().
#'
#' @return The same data.frame with two added columns:
#'   - Recurrent: "Yes" or "No"
#'   - n_samples: Number of supporting samples (if recurrent)
#' @export
#' @examples
#' input <- data.frame(
#'ID = c("S1", "S2", "S3", "S4"),
#'seqnames = c("chr1", "chr1", "chr1", "chr2"),
#'start = c(100, 120, 500, 100),
#'end = c(150, 170, 550, 150),
#'n_mendelian_error = c(10, 20, 5, 200)
#')
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
markRecurrentRegions <- function(df, recurrent_regions) {
  
  # --- Detect coordinate columns automatically ---
  if (all(c("start", "end") %in% colnames(df))) {
    start_col <- "start"
    end_col   <- "end"
  } else if (all(c("min_start", "max_end") %in% colnames(df))) {
    start_col <- "min_start"
    end_col   <- "max_end"
  } else {
    stop("Input must contain either (start, end) or (min_start, max_end) columns.")
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
    seqnames = df$seqnames,
    ranges = IRanges::IRanges(
      start = as.numeric(df[[start_col]]),
      end = as.numeric(df[[end_col]])
    ),
    ID = df$ID
  )
  
  
  overlaps <- GenomicRanges::findOverlaps(gr, recurrent_regions)
  
  if (length(overlaps) > 0) {
    df$Recurrent[S4Vectors::queryHits(overlaps)] <- "Yes"
    df$n_samples[S4Vectors::queryHits(overlaps)] <-
      recurrent_regions$n_samples[S4Vectors::subjectHits(overlaps)]
  }
  rownames(df) <- NULL
  return(df)
}
