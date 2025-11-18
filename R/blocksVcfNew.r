#' Collapse contiguous variants into genomic blocks based on HMM states
#'
#' This function creates block-level representations from a `CollapsedVCF`.
#' It performs the following steps:
#' 1. Identifies consecutive variants sharing the same inferred HMM state (`states`) and groups them into blocks.
#' 2. Extracts genomic coordinates, number of variants, HMM state, and genotype encodings (`geno_coded`) for each block.
#' 3. Optionally computes per-block depth ratios relative to a chromosome-wide mean (`total_mean`) from a specified FORMAT field (e.g., `DP` or `AD`) if provided.
#'
#' @param largeCollapsedVcf A `CollapsedVCF` containing `states` and `geno_coded` in `mcols()`.
#' @param add_ratios Logical; if `TRUE`, per-block ratios are computed and added to the output.
#' @param field_DP Optional character specifying the FORMAT depth field to use (e.g., `"DP"` or `"AD"`). Auto-selected if `NULL`.
#' @param total_mean Optional numeric vector containing chromosome-wide mean depth per sample, used to compute ratios.
#' @param ratio_cols Character vector of column names to assign to the ratio output when `total_mean` is provided.
#'
#' @return A `data.frame` with one row per block, containing chromosome, start, end, HMM state, number of SNPs,
#'         genotype encodings, and optionally per-block ratios relative to `total_mean`.
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
        ratio_matrix <- means / total_mean
        colnames(ratio_matrix) <- ratio_cols
        df <- cbind(df, as.data.frame(ratio_matrix))
      } 
    }
  }

  return(df)
}