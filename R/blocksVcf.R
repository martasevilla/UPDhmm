#' Collapse contiguous variants with the same state into blocks
#'
#' This function takes a data.table produced by `asDfVcf` and
#' identifies contiguous variants that share the same inferred state (`group`). 
#' These contiguous variants are then merged into "blocks", summarizing their start
#' and end positions, number of variants, genotype codings, and optionally 
#' summed quality/depth metrics.
#'
#' @param df data.table resulting from the `as_df_vcf` function.
#' @return data.table containing one row per block
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
    sum_matrix <- data.table::as.data.table(rowsum(quality_matrix, group = block_idx, na.rm = TRUE))
    
    # Count non-NA values per block
    count_matrix <- data.table::as.data.table(rowsum(1L * (!is.na(quality_matrix)), group = block_idx))
    
    # Rename columns for clarity
    colnames(sum_matrix)   <- paste0("total_sum_", quality_cols)
    colnames(count_matrix) <- paste0("total_count_", quality_cols)
  }

  ## --------------------------------------------------------------
  ## Build the resulting data.table
  ## --------------------------------------------------------------
  result <- data.table::data.table(
    ID         = sample_ID,
    seqnames   = seqnames_block,
    start      = start_block,
    end        = end_block,
    group      = r$values,
    n_snps     = r$lengths,
    geno_coded = geno_lists
  )

  ## Attach quality metrics if present
  if (!is.null(sum_matrix)) {
    result <- cbind(result, sum_matrix, count_matrix)
  }

  return(result)
}
