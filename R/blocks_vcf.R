#' Function to simplify contiguous variants with the same state into blocks.
#'
#' @param df Dataframe resulting from the `as_df_vcf` function.
#'
#' @return Dataframe containing information on the chromosome, start
#' position of the block, end position of the block, and predicted state.

blocks_vcf <- function(df) {
  setDT(df)
  df$n_snps_raw <- rleid(df$group)

  simplified_df <- df[, .(start = min(start), end = max(end), group = unique(group), seqnames = unique(seqnames), n_snps = .N), by = n_snps_raw]
  simplified_df<-simplified_df[,-c("n_snps_raw")]
  return(simplified_df)
}
