#' Function to transform a large collapsed VCF into a dataframe with
#' predicted states, including chromosome, start position, end position
#' and metadata.
#'
<<<<<<< HEAD
#' @param largeCollapsedVcf Name of the large collapsed VCF file.
#' @param genotypes Possible GT formats and its correspondence with the hmm
#'
#' @return dataframe
asDfVcf <- function(largeCollapsedVcf, genotypes) {
    genotypes_coded <- paste0(
        genotypes[VariantAnnotation::geno(largeCollapsedVcf)$GT[, "father"]],
        genotypes[VariantAnnotation::geno(largeCollapsedVcf)$GT[, "mother"]],
        genotypes[VariantAnnotation::geno(largeCollapsedVcf)$GT[, "proband"]]
    )
    vcf <- data.frame(
        start = GenomicRanges::start(largeCollapsedVcf),
        end = GenomicRanges::end(largeCollapsedVcf),
        group = S4Vectors::mcols(largeCollapsedVcf)$states,
        seqnames = as.character(GenomicRanges::seqnames(largeCollapsedVcf)),
        genotype = genotypes_coded
    )
=======
#' @param largecollapsedVcf Name of the large collapsed VCF file.
#' @param genotypes Possible GT formats and its correspondency with the hmm
#'
#' @return dataframe
asDfVcf <- function(largecollapsedVcf = NULL,genotypes = NULL) {
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
genotype = genotypes_coded)
>>>>>>> upstream/devel
}
