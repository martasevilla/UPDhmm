#' Function to transform a large collapsed VCF into a dataframe,
#' incorporating predicted states along with the log-likelihood
#' ratio and p-value.
#'
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



