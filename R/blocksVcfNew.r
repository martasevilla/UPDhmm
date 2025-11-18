#' Collapse contiguous variants with the same HMM state into blocks
#'
#' Creates genomic blocks directly from a `CollapsedVCF` object by grouping
#' consecutive variants that share the same inferred HMM state (`states`).
#' For each block, the function reports its genomic range, number of variants,
#' HMM state, and the list of genotype encodings (`geno_coded`). Optionally,
#' it also computes mean depth/quality metrics per sample.
#'
#' @param vcf A `CollapsedVCF` containing `states` and `geno_coded` in `mcols()`.
#' @param add_ratios Logical; if TRUE, compute mean DP/AD metrics per block.
#' @param field_DP Optional depth field to use (e.g. "DP", "AD"). Auto-detected if NULL.
#'
#' @return A data.frame with one row per block containing: chromosome, start, end,
#'   HMM state, number of SNPs, genotype encodings, and (optionally) mean depth/quality.
#'

blocksVcfNew <- function(largeCollapsedVcf, add_ratios = FALSE, field_DP = NULL, total_mean = NULL, ratio_cols = c("ratio_proband", "ratio_mother", "ratio_father")) {
  
  mcol <- S4Vectors::mcols(largeCollapsedVcf)
  coldata <- SummarizedExperiment::colData(largeCollapsedVcf)
  geno <- VariantAnnotation::geno(largeCollapsedVcf)

  states <- mcol$states
  geno_coded <- mcol$geno_coded
  seqnames_chr <- as.character(GenomicRanges::seqnames(largeCollapsedVcf))
  start_pos <- GenomicRanges::start(largeCollapsedVcf)
  end_pos <- GenomicRanges::end(largeCollapsedVcf)
  sample_ID <- coldata$ID[1]

  r <- rle(states)
  n_blocks <- length(r$lengths)
  ends_idx <- cumsum(r$lengths)
  starts_idx <- ends_idx - r$lengths + 1

  df <- data.frame(
    ID       = sample_ID,
    seqnames = seqnames_chr[starts_idx],
    start    = start_pos[starts_idx],
    end      = end_pos[ends_idx],
    group    = r$values,
    n_snps   = r$lengths,
    geno_coded = I(split(geno_coded, rep(seq_len(n_blocks), r$lengths)))
  )
  
  if (add_ratios) {
    dp_field <- if (!is.null(field_DP) && field_DP %in% names(geno)) { field_DP 
                } else if ("DP" %in% names(geno)) { "DP" 
                } else if ("AD" %in% names(geno)) { "AD" 
                } else { NULL }
    
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

      means <- sums / pmax(counts, 1)
      means[counts == 0] <- NA

      if (!is.null(total_mean)) {
        # Compute ratios block / total
        ratio_matrix <- means / total_mean
        colnames(ratio_matrix) <- ratio_cols
        df <- cbind(df, as.data.frame(ratio_matrix))
      } 
    }
  }

  return(df)
}