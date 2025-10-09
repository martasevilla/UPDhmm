#' Function to transform a large collapsed VCF into a dataframe,
#' incorporating predicted states along with the log-likelihood
#' ratio and p-value.
#'
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
  #Create a granges with the interval
  gRanges_block <- GenomicRanges::makeGRangesFromDataFrame(
    filtered_def_blocks_states,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = c("seqnames"),
    start.field = "start",
    end.field = c("end")
  )
  #Overlap it with the raw largeCollapsedVcf to extract the genotypes
  overlaps <- IRanges::subsetByOverlaps(largeCollapsedVcf, gRanges_block)
  #Transform the gentoypes
  genotypes_uncoded <- VariantAnnotation::geno(overlaps)$GT
  genotypes_coded <- c(
    paste0(
      genotypes[genotypes_uncoded[, "father"]],
      genotypes[genotypes_uncoded[, "mother"]],
      genotypes[genotypes_uncoded[, "proband"]]
    )
  )
  
  
  
  #Create a function to measure the cumulative probability of having that 
  #state at given genotype in that chain of observations
  prob_state_chain <- function(state, observed_sequence, hmm_model) {
    N <- length(hmm_model$States)
    Tn <- length(observed_sequence)
    
    # Find the index of the state
    state_index <- match(state, hmm_model$States)
    
    # Initialization
    alpha <- numeric(Tn)
    alpha[1] <- hmm_model$startProbs[state_index] * 
      hmm_model$emissionProbs[state_index, observed_sequence[1]]
    
    # Recursion
    for (t in 2:Tn) {
      alpha[t] <- hmm_model$transProbs[state_index, state_index] * 
        hmm_model$emissionProbs[state_index, observed_sequence[t]]
    }
    
    # Probability of the observed sequence given the state
    return(alpha)
  }
  
  #Apply to a normal state
  prob_matrix_normal <- prob_state_chain("normal", 
                                         genotypes_coded, 
                                         hmm_model = hmm)
  #Apply to the predicted state
  prob_matrix_other <- prob_state_chain(S4Vectors::mcols(gRanges_block)$group, 
                                        genotypes_coded, hmm_model = hmm)
  
  #Calculate the overal probability of the whole interval for normal state and 
  #predicted state
  log_likelihood_normal <- sum(log(prob_matrix_normal))
  log_likelihood_other <- sum(log(prob_matrix_other))
  
  #Calculate the likelihood ratio between the predicted state and normal state
  #and its related p-value
  overall_log_odds_ratio <- -2 * (log_likelihood_normal - log_likelihood_other)
  p_value <- stats::pchisq(overall_log_odds_ratio, 1, lower.tail = FALSE)
  
  #Assign it into metadata columns
  S4Vectors::mcols(gRanges_block)$log_likelihood <- overall_log_odds_ratio
  S4Vectors::mcols(gRanges_block)$p_value <- p_value
  
  #Now, retrieve the count of Mendelian errors (genotypic combinations that 
  #are inconsistent with Mendelian inheritance principles).
  emission_probs <- hmm$emissionProbs["normal",]
  mendelian_error_values<-names(emission_probs[emission_probs == 0.00001])
  n_mendelian_error<-sum(genotypes_coded %in% mendelian_error_values)
  
  #Add it as a metadata column
  S4Vectors::mcols(gRanges_block)$n_mendelian_error <- n_mendelian_error
  #create a dataframe with the log_likelihood, p-value and n_mendelian errors
  df <- as.data.frame(gRanges_block)
  df <- df[, !(names(df) %in% c("strand", "width"))]
  df <- df[, c("ID", setdiff(names(df), "ID"))]
  return(df)
}
