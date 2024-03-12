#' Function to transform a large collapsed VCF into a dataframe with
#' predicted states, including chromosome, start position, end position
#' and metadata.
#'
#' @param largecollapsedVcf Name of the large collapsed VCF file.
#' @param genotypes Possible GT formats and its correspondency with the hmm
#'
#' @return dataframe
as_df_vcf <- function(largecollapsedVcf = NULL,genotypes = NULL) {


  genotypes_coded <- paste0(
  genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "father"]],
  genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "mother"]],
  genotypes[VariantAnnotation::geno(largecollapsedVcf)$GT[, "proband"]]
  )
  vcf <- data.frame(
  start = GenomicRanges::start(largecollapsedVcf),
  end = GenomicRanges::end(largecollapsedVcf),
  group = S4Vectors::mcols(largecollapsedVcf)$states,
  seqnames = as.character(GenomicRanges::seqnames(largecollapsedVcf)),
  genotype = genotypes_coded

  )
}
