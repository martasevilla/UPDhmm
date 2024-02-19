#' Function for simplify into blocks contiguous variants with same state
#'
#' @param df df result from as_df_vcf function
#'
#' @return dataframe
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
