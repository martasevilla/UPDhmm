#' Collapse contiguous variants into genomic blocks based on HMM states
#'
#' This internal helper function creates block-level representations from a 
#' \code{CollapsedVCF} object. It identifies consecutive variants sharing the 
#' same inferred HMM state (\code{states}) and groups them into contiguous blocks.
#' For each block, it summarizes genomic coordinates, number of variants, HMM state, 
#' genotype encodings (\code{geno_coded}), and optionally computes per-block depth 
#' ratios relative to a chromosome-wide mean.
#'
#' @param largeCollapsedVcf A \code{CollapsedVCF} object containing \code{states} 
#'   and \code{geno_coded} in \code{mcols()}.
#' @param add_ratios Logical; default \code{FALSE}. If \code{TRUE}, per-block ratios
#'   are computed from a DP-like field.
#' @param field_DP Optional character specifying the FORMAT field to use for depth 
#'   (e.g., \code{"DP"} or \code{"AD"}). If \code{NULL}, the function auto-selects.
#' @param total_mean Optional numeric vector with chromosome-wide mean depth per 
#'   sample, used to compute block-level ratios.
#' @param ratio_cols Character vector of column names to assign to the ratio output 
#'   when \code{total_mean} is provided. Default: \code{c("ratio_proband", "ratio_mother", "ratio_father")}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Identifies consecutive variants with the same \code{states} using run-length encoding.
#'   \item Constructs a block-level \code{data.frame} containing:
#'         \itemize{
#'           \item \code{seqnames}, \code{start}, \code{end} – genomic coordinates
#'           \item \code{group} – HMM state of the block
#'           \item \code{n_snps} – number of variants in the block
#'           \item \code{geno_coded} – list of numeric genotype codes for the block
#'         }
#'   \item Optionally computes per-block depth ratios relative to \code{total_mean} 
#'         if \code{add_ratios = TRUE} and a valid depth field is present.
#' }
#'
#' @return A \code{data.frame} with one row per block, containing:
#'   \itemize{
#'     \item \code{ID} – sample identifier
#'     \item \code{seqnames}, \code{start}, \code{end} – genomic coordinates
#'     \item \code{group} – HMM state of the block
#'     \item \code{n_snps} – number of variants in the block
#'     \item \code{geno_coded} – list of numeric genotype codes per block
#'     \item Optional ratio columns relative to \code{total_mean}
#'   }
#'
blocksVcfNew <- function(largeCollapsedVcf, add_ratios = FALSE, field_DP = NULL, total_mean = NULL, ratio_cols = c("ratio_proband", "ratio_mother", "ratio_father")) {
  
  # Extract metadata and genotype matrices
  mcol <- S4Vectors::mcols(largeCollapsedVcf)
  coldata <- SummarizedExperiment::colData(largeCollapsedVcf)
  geno <- VariantAnnotation::geno(largeCollapsedVcf)

  # Extract per-variant information for block construction
  states <- mcol$states
  geno_coded <- mcol$geno_coded
  seqnames_chr <- as.character(GenomicRanges::seqnames(largeCollapsedVcf))
  start_pos <- GenomicRanges::start(largeCollapsedVcf)
  end_pos <- GenomicRanges::end(largeCollapsedVcf)
  sample_ID <- coldata$ID[1]

  # Identify contiguous blocks of variants sharing the same HMM state
  r <- rle(states)
  n_blocks <- length(r$lengths)
  ends_idx <- cumsum(r$lengths)
  starts_idx <- ends_idx - r$lengths + 1

  # Construct initial data.frame with block information
  df <- data.frame(
    ID       = sample_ID,
    seqnames = seqnames_chr[starts_idx],
    start    = start_pos[starts_idx],
    end      = end_pos[ends_idx],
    group    = r$values,
    n_snps   = r$lengths,
    geno_coded = I(split(geno_coded, rep(seq_len(n_blocks), r$lengths)))
  )
  
  # Optionally compute per-block depth ratios
  if (add_ratios) {
    dp_field <- if (!is.null(field_DP) && field_DP %in% names(geno)) { field_DP } 
                else if ("DP" %in% names(geno)) { "DP" } 
                else if ("AD" %in% names(geno)) { "AD" } 
                else { NULL }
    
    if (!is.null(dp_field)) {
      if (dp_field == "AD") {
        quality_matrix <- apply(geno$AD, 2, function(col) {
          vapply(col, function(x) {if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)}, numeric(1))
        })
        
      } else {
        quality_matrix <- as.matrix(geno[[dp_field]])
      }
      
      block_idx <- rep(seq_len(n_blocks), r$lengths)

      # Sum quality values per block, ignoring NA
      sums <- rowsum(quality_matrix, group = block_idx, na.rm = TRUE)
    
      # Count non-NA values per block
      counts <- rowsum(1L * (!is.na(quality_matrix)), group = block_idx)

      # Compute means per block
      means <- sums / pmax(counts, 1)
      means[counts == 0] <- NA
      
      if (!is.null(total_mean)) {
        # Compute ratios relative to total_mean
        ratio_matrix <- sweep(means, 2, total_mean, FUN = "/")
        colnames(ratio_matrix) <- ratio_cols
        df <- cbind(df, as.data.frame(ratio_matrix))
      } 
    }
  }

  return(df)
}