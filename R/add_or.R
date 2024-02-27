#' Function to transform a large collapsed VCF into a dataframe, incorporating
#' predicted states along with the log-likelihood ratio and p-value.
#'
#' @param split_vcf_df Input VCF file in dataframe format.
#' @param start_coord Start chromosomal position for each block.
#' @param end_coord End chromosomal position for each block.
#' @param seqnames Name of the chromosome for the blocks (e.g., chr1 or 1).
#' @param group Hidden state assigned by UPDHmm.
#' @return Dataframe containing the transformed information.

add_or <- function(start_coord, end_coord, seqnames, group, split_vcf_df = NULL) {
  if (is.null(split_vcf_df)) {
    split_vcf_df <- split_vcf_df
  }

  # Extract relevant genotypes from vcf_df
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
  #Extract p-value
  p_value <- stats::pchisq(overall_log_odds_ratio, 1, lower.tail = FALSE)

  # Return a vector with log_OR and p_value
  result <- c(log_OR = overall_log_odds_ratio, p_value = p_value)
  names(result) <- c("log_OR", "p_value")

  return(result)
}
