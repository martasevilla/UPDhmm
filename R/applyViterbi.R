#' Apply the hidden Markov model using the Viterbi algorithm.
#'
#'
#' @param largeCollapsedVcf input vcf file
#' @param hmm Hidden Markov Model used to infer the events
#' @param genotypes Possible GT formats and its correspondence with the hmm
#' @return largeCollapsedVcf

applyViterbi <- function(largeCollapsedVcf, hmm) {

  ## Retrieve the precomputed numeric genotypes (observations) from the VCF metadata
  geno_coded_values <- S4Vectors::mcols(largeCollapsedVcf)$geno_coded

  ## Run the Viterbi algorithm to infer the most likely sequence of hidden states
  S4Vectors::mcols(largeCollapsedVcf)$states <- HMM::viterbi(hmm, geno_coded_values )

  ## Return the updated VCF object with inferred states stored in metadata
  return(largeCollapsedVcf)
}
