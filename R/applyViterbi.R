#' Apply the Hidden Markov Model using the Viterbi algorithm
#'
#' This internal helper function takes a \code{CollapsedVCF} object and a 
#' pre-defined Hidden Markov Model (HMM), applies the Viterbi algorithm to infer 
#' the most likely sequence of hidden states for each variant, and stores the 
#' resulting state sequence inside the VCF metadata (\code{mcols}).
#'
#' The function assumes that a pre-computed numeric vector of genotype codes 
#' (\code{geno_coded}) is already present in the metadata columns of the VCF, 
#' typically generated during earlier preprocessing steps.
#'
#' @param largeCollapsedVcf A \code{CollapsedVCF} object containing the variants 
#'   for a chromosome.
#' @param hmm A Hidden Markov Model (HMM) object used to infer hidden states 
#'   through the Viterbi algorithm.
#'
#' @return The input \code{CollapsedVCF} object with an additional metadata 
#'   column \code{states} containing the inferred hidden state at each variant.
#'
applyViterbi <- function(largeCollapsedVcf, hmm) {

  ## Retrieve the precomputed numeric genotypes (observations) from the VCF metadata
  geno_coded_values <- S4Vectors::mcols(largeCollapsedVcf)$geno_coded

  ## Run the Viterbi algorithm to infer the most likely sequence of hidden states
  S4Vectors::mcols(largeCollapsedVcf)$states <- HMM::viterbi(hmm, geno_coded_values )

  ## Return the updated VCF object with inferred states stored in metadata
  return(largeCollapsedVcf)
}
