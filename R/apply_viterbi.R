#' Apply the hidden Markov model using the Viterbi algorithm.
#'
#'
#' @param largecollapsedVcf input vcf file
#' @param hmm Hidden Markov Model used to infer the events
#' @param genotypes Possible GT formats and its correspondency with the hmm
#' @return largecollapsedVcf

apply_viterbi<- function(largecollapsedVcf=NULL,hmm=NULL,genotypes=NULL) {

# First, assign the names of the samples according to the trio.
# It is crucial to maintain the order for applying the Viterbi algorithm.

  vector_samples <- colnames(largecollapsedVcf)

# Now, transform genotypes into numerical codes and apply Viterbi algorithm.
  # The results will be stored as a metadata column in the object.
  genotypes_coded<<-c(paste0(
    genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "father"]],
    genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "mother"]],
    genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "proband"]]
  ))

  states<-HMM::viterbi(hmm,genotypes_coded)
  S4Vectors::mcols(largecollapsedVcf)$states <-states
  return(largecollapsedVcf)

}
