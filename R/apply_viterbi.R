#' Apply the hidden Markov model using the Viterbi algorithm.
#'
#'
#' @param largecollapsedVcf input vcf file
#'
#' @return largecollapsedVcf

apply_viterbi <- function(largecollapsedVcf) {
  # First, assign the names of the samples according to the trio.
  # It is crucial to maintain the order for applying the Viterbi algorithm.
  vector_samples <- rownames(SummarizedExperiment::colData(largecollapsedVcf))
  # Now, transform genotypes into numerical codes and apply the Viterbi algorithm.
  # The results will be stored as a metadata column in the object.
  genotypes <- c(
    "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
    "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
  )
  utils::data("hmm")
  hmm <- hmm
  GenomicRanges::elementMetadata(largecollapsedVcf)$metadata <-
    HMM::viterbi(
      hmm,
      c(paste0(
        genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "father"]],
        genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "mother"]],
        genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "proband"]]
      ))
    )
  return(largecollapsedVcf)
}
