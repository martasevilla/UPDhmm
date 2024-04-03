#' Apply the hidden Markov model using the Viterbi algorithm.
#'
#'
<<<<<<< HEAD
#' @param largeCollapsedVcf input vcf file
#' @param hmm Hidden Markov Model used to infer the events
#' @param genotypes Possible GT formats and its correspondence with the hmm
#' @return largeCollapsedVcf

applyViterbi <-
  function(largeCollapsedVcf, hmm, genotypes) {
    ## First, assign the names of the samples according to the trio.
    ## It is crucial to maintain the order for applying the Viterbi algorithm.
        vector_samples <- colnames(largeCollapsedVcf)

    ## Now, transform genotypes into numerical codes and apply Viterbi algorithm.
    ## The results will be stored as a metadata column in the object.
        genotypes_coded <- c(paste0(
          genotypes[VariantAnnotation::geno(largeCollapsedVcf)$GT[, "father"]],
          genotypes[VariantAnnotation::geno(largeCollapsedVcf)$GT[, "mother"]],
          genotypes[VariantAnnotation::geno(largeCollapsedVcf)$GT[, "proband"]]
        ))

        states <- HMM::viterbi(hmm, genotypes_coded)
        S4Vectors::mcols(largeCollapsedVcf)$states <- states
        return(largeCollapsedVcf)
    }
=======
#' @param largecollapsedVcf input vcf file
#' @param hmm Hidden Markov Model used to infer the events
#' @param genotypes Possible GT formats and its correspondency with the hmm
#' @return largecollapsedVcf

applyViterbi <-
function(largecollapsedVcf = NULL,hmm = NULL,genotypes = NULL)
{
## First, assign the names of the samples according to the trio.
## It is crucial to maintain the order for applying the Viterbi algorithm.
vector_samples <- colnames(largecollapsedVcf)

## Now, transform genotypes into numerical codes and apply Viterbi algorithm.
## The results will be stored as a metadata column in the object.
genotypes_coded <- c(paste0(
genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "father"]],
genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "mother"]],
genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "proband"]]
))

states <- HMM::viterbi(hmm,genotypes_coded)
S4Vectors::mcols(largecollapsedVcf)$states <- states
return(largecollapsedVcf)
}
>>>>>>> upstream/devel
