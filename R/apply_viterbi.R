#' This function is for applying the hidden-markov-model using the
#' viterbi algorithm
#' IMPORTANT: vector_samples should be previously created
#'
#' @param largecollapsedVcf input vcf file
#'
#' @return largecollapsedVcf

apply_viterbi <- function(largecollapsedVcf) {
  # first assign the names of the samples, according to the trio,
  # it is important the order for appliying Viterbi algortihm
  vector_samples <- rownames(SummarizedExperiment::colData(largecollapsedVcf))
  # now, tranfsofr genotypes into numerical code, and apply viterbi alogortihm.
  # The result would be stored as a metatadata column in the object
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
