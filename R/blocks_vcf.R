#' Function to simplify contiguous variants with the same state into
#' blocks.
#'
#' @param df Dataframe resulting from the `as_df_vcf` function.
#' @return Dataframe containing information on the chromosome, start
#' #' position of the block, end position of the block, and predicted state.
#' @importFrom magrittr "%>%"
#' @import dplyr


blocks_vcf <- function(df) {
  n_snps_raw<-NULL
  end<-NULL
  group<-NULL
  seqnames<-NULL
  start<-NULL
  df <- df %>%
  dplyr::mutate(n_snps_raw = cumsum(c(TRUE,
                                      group[-1L] != group[-length(group)])))

  simplified_df <- df %>%
  dplyr::group_by(n_snps_raw) %>%
  dplyr::summarise(
    start = min(start),
    end = max(end),
    group = unique(group),
    seqnames = unique(seqnames),
    n_snps = n()
  )


  simplified_df$n_snps_raw<-NULL

  return(simplified_df)
}
