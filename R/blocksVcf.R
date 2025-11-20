#' Collapse contiguous variants with the same inferred state into blocks
#'
#' This function identifies consecutive variants in a data.frame
#' (produced by \code{asDfVcf()}) that share the same inferred state 
#' and merges them into contiguous "blocks". Each block summarizes its start 
#' and end positions, the number of variants it contains, genotype codings, and 
#' optionally summed quality or depth metrics if the columns exist in the input.
#'
#' @param df A data.frame produced by \code{asDfVcf()}.
#'
#' @return A data.frame with one row per block. Columns include:
#' \itemize{
#'   \item ID – sample identifier
#'   \item seqnames, start, end – genomic coordinates
#'   \item group – inferred state of the block
#'   \item n_snps – number of variants in the block
#'   \item geno_coded – list of numeric genotype codes per block
#'   \item summed and counted quality/depth metrics per block (if present)
#' }
#'
blocksVcf <- function(df) {

  ## Identify contiguous blocks of the same state using run-length encoding
  r <- rle(df$group)
  n_blocks <- length(r$lengths) # number of contiguous blocks
  ends_idx <- cumsum(r$lengths) # ending index of each block in the original df
  starts_idx <- ends_idx - r$lengths + 1L # starting index of each block
  block_idx <- rep(seq_len(n_blocks), r$lengths) # vector mapping each variant to a block

  sample_ID <- unique(df$ID)[1]

  # Extract block-level genomic information
  seqnames_block <- df$seqnames[starts_idx]  # chromosome per block
  start_block    <- df$start[starts_idx]     # start position of block
  end_block      <- df$end[ends_idx]         # end position of block

  ## Collect genotype codings for each block
  geno_lists <- split(df$geno_coded, block_idx)

  ## --------------------------------------------------------------
  ## Optional: sum and count quality metrics per block
  ## --------------------------------------------------------------
  quality_cols <- intersect(c("quality_proband", "quality_mother", "quality_father"), names(df))
  sum_matrix <- count_matrix <- NULL

  if (length(quality_cols) > 0L) {

    quality_matrix <- as.matrix(df[, quality_cols, drop = FALSE])

    # Sum quality values per block, ignoring NA
    sum_matrix <- as.data.frame(rowsum(quality_matrix, group = block_idx, na.rm = TRUE))
    
    # Count non-NA values per block
    count_matrix <- as.data.frame(rowsum(1L * (!is.na(quality_matrix)), group = block_idx))
    
    # Rename columns for clarity
    colnames(sum_matrix)   <- paste0("total_sum_", quality_cols)
    colnames(count_matrix) <- paste0("total_count_", quality_cols)

  }

  ## --------------------------------------------------------------
  ## Build the resulting data.frame
  ## --------------------------------------------------------------
  result <- data.frame(
    ID         = sample_ID,
    seqnames   = seqnames_block,
    start      = start_block,
    end        = end_block,
    group      = r$values,
    n_snps     = r$lengths,
    geno_coded = I(geno_lists)
  )
  
  ## Attach quality metrics if present
  if (!is.null(sum_matrix)) {
    result <- cbind(result, sum_matrix, count_matrix)
  }

  return(result)
}
