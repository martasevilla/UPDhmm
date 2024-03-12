#' Function to simplify contiguous variants with the same state into
#' blocks.
#'
#' @param df data.frame resulting from the `as_df_vcf` function.
#' @return data.frame containing information on the chromosome, start
#' #' position of the block, end position of the block, and predicted state.


blocks_vcf <- function(df) {

  # Add a column for n_snps_raw
  df$n_snps_raw <- base::cumsum(c(TRUE, df$group[-1L] != df$group[-length(df$group)]))

  # Create an empty data frame to store the aggregated results
  simplified_df <- data.frame(start = integer(),
                              end = integer(),
                              group = character(),
                              seqnames = character(),
                              n_snps = integer(),
                              stringsAsFactors = FALSE)

  # Loop through unique n_snps_raw values
  for (n in cumsum::unique(df$n_snps_raw)) {
    subset_df <- df[df$n_snps_raw == n, ]
    new_row <- data.frame(
      start = min(subset_df$start),
      end = max(subset_df$end),
      group = unique(subset_df$group),
      seqnames = unique(subset_df$seqnames),
      n_snps = nrow(subset_df),
      stringsAsFactors = FALSE
    )
    simplified_df <- rbind(simplified_df, new_row)
  }

  # Remove the n_snps_raw column
  simplified_df$n_snps_raw <- NULL
  return(simplified_df)
}
