#' Function for transform a largecollapsedVcf into a dataframe with
#' predicted states (only with chr,start,end and metadatacolumn)
#'
#' @param split_vcf_df input vcf file
#' @param start_coord start chromosomal position of every block
#' @param end_coord  end chromosomal position of every block
#' @param seqnames name of the chromosome of the blocks (e.g. chr1 or 1)
#' @param group hidden state assigned by  UPDHmm
#' @return dataframe

add_or <- function(start_coord, end_coord, seqnames, group, split_vcf_df = NULL) {
  if (is.null(split_vcf_df)) {
    split_vcf_df <- split_vcf_df
  }

  # Extract relevant genotypes from as_vcf_df
  genotypes <-
    split_vcf_df[[seqnames]]$genotype[split_vcf_df[[seqnames]]$start >= start_coord &
      split_vcf_df[[seqnames]]$end <= end_coord]

  # Run forward algorithm for the block
  utils::data("hmm")
  hmm <- hmm
  forward_matrix <- HMM::forward(hmm, genotypes)

  log_likelihood_normal <- utils::tail(forward_matrix["normal", ], 1)
  log_likelihood_other <- utils::tail(forward_matrix[group, ], 1)
  overall_log_odds_ratio <- -2 * (log_likelihood_normal - log_likelihood_other)
  p_value <- stats::pchisq(overall_log_odds_ratio, 1, lower.tail = FALSE)

  # Return a named vector with log_OR and p_value
  result <- c(log_OR = overall_log_odds_ratio, p_value = p_value)
  names(result) <- c("log_OR", "p_value")

  return(result)
}
