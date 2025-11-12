#' Apply the hidden Markov model using the Viterbi algorithm.
#'
#' This function takes a collapsed VCF object and a pre-defined Hidden Markov Model (HMM),
#' applies the Viterbi algorithm to infer the most likely sequence of hidden states for 
#' each sample, and stores the resulting states as a new metadata column in the VCF object.
#'
#' @param largeCollapsedVcf input vcf file
#' @param hmm Hidden Markov Model used to infer the events
#' @return largeCollapsedVcf
#'

applyViterbi <- function(largeCollapsedVcf, hmm) {

  ## Retrieve the precomputed numeric genotypes (observations) from the VCF metadata
  geno_matrix <- S4Vectors::mcols(largeCollapsedVcf)$geno_coded

  ## Run the Viterbi algorithm to infer the most likely sequence of hidden states
  S4Vectors::mcols(largeCollapsedVcf)$states <- HMM::viterbi(hmm, geno_matrix)

  ## Return the updated VCF object with inferred states stored in metadata
  return(largeCollapsedVcf)
}
