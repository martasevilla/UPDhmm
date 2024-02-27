#' Function to simplify contiguous variants with the same state into blocks.
#'
#' @param df Dataframe resulting from the `as_df_vcf` function.
#'
#' @return Dataframe containing information on the chromosome, start position of the block, end position of the block, and predicted state.

blocks_vcf <- function(df) {
  df$n_snps <- 1

  for (i in 2:nrow(df)) {
    if (df$group[i] == df$group[i - 1]) {
      df$n_snps[i] <- df$n_snps[i - 1] + 1
    }
  }

  simplified_df <- data.frame(
    start = tapply(df$start, cumsum(c(TRUE, diff(df$n_snps) != 1)), min),
    end = tapply(df$end, cumsum(c(TRUE, diff(df$n_snps) != 1)), max),
    group = tapply(
      as.character(df$group), cumsum(c(
        TRUE,
        diff(df$n_snps) != 1
      )),
      unique
    ),
    n_snps = tapply(df$n_snps, cumsum(c(TRUE, diff(df$n_snps) != 1)), max)
  )

  simplified_df$seqnames <- rep(df$seqnames[1], nrow(simplified_df))
  return(simplified_df)
}
