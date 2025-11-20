#' Apply the Hidden Markov Model using the Viterbi algorithm
#'
#' This function takes a VCF object and a pre-defined Hidden Markov Model (HMM), 
#' applies the Viterbi algorithm to infer the most likely sequence of hidden states 
#' for each variant, and stores the resulting state sequence inside the VCF metadata.
#'
#' It expects that the metadata column *geno_coded* exists, containing a numeric
#' encoding of the genotypes. This column is automatically generated when the VCF
#' is processed with \code{vcfCheck()} from the UPDhmm package.
#'
#' @param largeCollapsedVcf A CollapsedVCF object. 
#' @param hmm A Hidden Markov Model object.
#'
#' @return The input CollapsedVCF object updated with a new metadata column *states*,
#'   which contains the inferred hidden state for each variant.
#'
applyViterbi <- function(largeCollapsedVcf, hmm) {

  ## Retrieve the precomputed numeric genotypes (observations) from the VCF metadata
  geno_coded_values <- S4Vectors::mcols(largeCollapsedVcf)$geno_coded

  ## Run the Viterbi algorithm to infer the most likely sequence of hidden states
  S4Vectors::mcols(largeCollapsedVcf)$states <- HMM::viterbi(hmm, geno_coded_values )

  ## Return the updated VCF object with inferred states stored in metadata
  return(largeCollapsedVcf)
}
