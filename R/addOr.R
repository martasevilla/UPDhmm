#' Function to transform a large collapsed VCF into a dataframe,
#' incorporating predicted states along with the log-likelihood
#' ratio and p-value.
#'
<<<<<<< HEAD
#' @param largeCollapsedVcf Input VCF file
#' @param filtered_def_blocks_states data.frame object containing the blocks
#' @param hmm Hidden Markov Model used to infer the events. The format should
#' adhere to the general HMM format from HMM package with a series of elements:
#' 1. The hidden states names in the "States" vector.
#' 2. All possible observations in the "Symbols" vector.
#' 3. Start probabilities of every hidden state in the "startProbs" vector.
#' 4. Transition probabilities matrix of the hidden states in "transProbs".
#' 5. Probabilities associated between every hidden state and all possible
#' observations in the "emissionProbs" matrix.
#' @param genotypes Possible GT formats and its correspondence with the hmm
#' @return data.frame containing the transformed information.

addOr <- function(
    filtered_def_blocks_states,
    largeCollapsedVcf,
    hmm,
    genotypes
) {
  gRanges_block <- GenomicRanges::makeGRangesFromDataFrame(
    filtered_def_blocks_states,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = c("seqnames"),
    start.field = "start",
    end.field = c("end")
  )
  
  overlaps <- IRanges::subsetByOverlaps(largeCollapsedVcf, gRanges_block)
  genotypes_uncoded <- VariantAnnotation::geno(overlaps)$GT
  genotypes_coded <- c(
    paste0(
      genotypes[genotypes_uncoded[, "father"]],
      genotypes[genotypes_uncoded[, "mother"]],
      genotypes[genotypes_uncoded[, "proband"]]
    )
  )
  
  
  
  
  forward_algorithm_state <- function(state, observed_sequence, hmm_model) {
    N <- length(hmm_model$States)
    Tn <- length(observed_sequence)
    
    # Find the index of the state
    state_index <- match(state, hmm_model$States)
    
    # Initialization
    alpha <- numeric(Tn)
    alpha[1] <- hmm_model$startProbs[state_index] * hmm_model$emissionProbs[state_index, observed_sequence[1]]
    
    # Recursion
    for (t in 2:Tn) {
      alpha[t] <- hmm_model$transProbs[state_index, state_index] * hmm_model$emissionProbs[state_index, observed_sequence[t]]
    }
    
    # Probability of the observed sequence given the state
    return(alpha)
  }
  
  forward_matrix_normal <- forward_algorithm_state("normal", genotypes_coded, hmm_model = hmm)
  forward_matrix_other <- forward_algorithm_state(S4Vectors::mcols(gRanges_block)$group, genotypes_coded, hmm_model = hmm)
  
  log_likelihood_normal <- sum(log(forward_matrix_normal))
  log_likelihood_other <- sum(log(forward_matrix_other))
  
  overall_log_odds_ratio <- -2 * (log_likelihood_normal - log_likelihood_other)
  p_value <- stats::pchisq(overall_log_odds_ratio, 1, lower.tail = FALSE)
  
  S4Vectors::mcols(gRanges_block)$log_OR <- overall_log_odds_ratio
  S4Vectors::mcols(gRanges_block)$p_value <- p_value
  
  df <- as.data.frame(gRanges_block)
  df <- df[, !(names(df) %in% c("strand", "width"))]
  return(df)
}
=======
#' @param largecollapsedVcf Input VCF file
#' @param filtered_def_blocks_states data.frame object containing the blocks
#' @param hmm Hidden Markov Model used to infer the events
#' @param genotypes Possible GT formats and its correspondency with the hmm
#' @return data.frame containing the transformed information.

addOr <-
function(filtered_def_blocks_states = NULL,
largecollapsedVcf = NULL, hmm = NULL, genotypes = NULL) {

gRanges_block <-
GenomicRanges::makeGRangesFromDataFrame(filtered_def_blocks_states,
keep.extra.columns = TRUE,
ignore.strand = TRUE,
seqnames.field = c("seqnames"),
start.field = "start",
end.field = c("end"))

overlaps <- IRanges::subsetByOverlaps(largecollapsedVcf, gRanges_block)
genotypes_uncoded <- VariantAnnotation::geno(overlaps)$GT
genotypes_coded <- c(
paste0(genotypes[genotypes_uncoded[, "father"]],
genotypes[genotypes_uncoded[, "mother"]],
genotypes[genotypes_uncoded[, "proband"]]))

forward_matrix <- HMM::forward(hmm, genotypes_coded)

log_likelihood_normal <-
utils::tail(forward_matrix["normal", ], 1)

log_likelihood_other <-
utils::tail(forward_matrix[S4Vectors::mcols(gRanges_block)$group, ], 1)

overall_log_odds_ratio <- -2 * (log_likelihood_normal - log_likelihood_other)
p_value <- stats::pchisq(overall_log_odds_ratio, 1, lower.tail = FALSE)

S4Vectors::mcols(gRanges_block)$log_OR <- overall_log_odds_ratio
S4Vectors::mcols(gRanges_block)$p_value <- p_value

df <- as.data.frame(gRanges_block)
df <- df[, !(names(df) %in% c("strand","width"))]
return(df)
}



>>>>>>> upstream/devel
