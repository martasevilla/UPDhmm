#' Function for transform a largecollapsedVcf into a dataframe with predicted
#' states (only with chr,start,end and metadatacolumn)
#'
#' @param largecollapsedVcf Vcf file
#'
#' @return dataframe
as_df_vcf <- function(largecollapsedVcf) {
  genotypes <- c(
    "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
    "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
  )
  vcf <- data.frame(
    start = GenomicRanges::start(largecollapsedVcf),
    end = GenomicRanges::end(largecollapsedVcf),
    group = GenomicRanges::elementMetadata(largecollapsedVcf)$metadata,
    seqnames = as.character(GenomicRanges::seqnames(largecollapsedVcf)),
    genotype =
      paste0(
        genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "father"]],
        genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "mother"]],
        genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "proband"]]
      )
  )
}
