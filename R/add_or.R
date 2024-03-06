#' Function to transform a large collapsed VCF into a dataframe, incorporating
#' predicted states along with the log-likelihood ratio and p-value.
#'
#' @param largecollapdsedVcf Input VCF file
#' @param filtered_def_blocks_states data.frame object containing the blocks
#' @param hmm Hidden Markov Model used to infer the events
#' @param genotypes Possible GT formats and its correspondency with the hmm
#' @return data.frame containing the transformed information.



add_or <- function(filtered_def_blocks_states=NULL,
                   largecollapsedVcf=NULL, hmm=NULL, genotypes=NULL) {

  gRanges_block<-
    GenomicRanges::makeGRangesFromDataFrame(filtered_def_blocks_states,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=TRUE,
                                            seqnames.field=c("seqnames"),
                                            start.field="start",
                                            end.field=c("end"))

  overlaps <- IRanges::subsetByOverlaps(largecollapsedVcf, gRanges_block)

  genotypes_uncoded <- VariantAnnotation::geno(overlaps)$GT

  genotypes_coded <- c(
    paste0(genotypes[genotypes_uncoded[, "father"]],
           genotypes[genotypes_uncoded[, "mother"]],
           genotypes[genotypes_uncoded[, "proband"]])
  )

  forward_matrix <- HMM::forward(hmm, genotypes_coded)

  log_likelihood_normal <-
    tail(forward_matrix["normal", ], 1)

  log_likelihood_other <-
    tail(forward_matrix[gRanges_block@elementMetadata$group, ], 1)


  overall_log_odds_ratio <- -2 * (log_likelihood_normal - log_likelihood_other)

  p_value <- pchisq(overall_log_odds_ratio, 1, lower.tail = FALSE)

  gRanges_block@elementMetadata$log_OR <- overall_log_odds_ratio
  gRanges_block@elementMetadata$p_value <- p_value


  df<-as.data.frame(gRanges_block)
  df <- df[, !(names(df) %in% "strand")]
  return(df)
}








