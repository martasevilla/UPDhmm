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

  ## Extract the genotype matrix from the vcf object
  geno_matrix <- S4Vectors::mcols(largeCollapsedVcf)$geno_coded

  ## Apply the Viterbi algorithm to infer the most likely hidden states and 
  S4Vectors::mcols(largeCollapsedVcf)$states <- HMM::viterbi(hmm, geno_matrix)

  ## Return the updated VCF object with the inferred states
  return(largeCollapsedVcf)
}
