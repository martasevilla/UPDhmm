#' Function to transform a large collapsed VCF into a dataframe with
#' predicted states, including chromosome, start position, end position
#' and metadata.
#'
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
        ID =  SummarizedExperiment::colData(largeCollapsedVcf)$ID[SummarizedExperiment::colData(largeCollapsedVcf)$Samples == 1],
        start = GenomicRanges::start(largeCollapsedVcf),
        end = GenomicRanges::end(largeCollapsedVcf),
        group = S4Vectors::mcols(largeCollapsedVcf)$states,
        seqnames = as.character(GenomicRanges::seqnames(largeCollapsedVcf)),
        genotype = genotypes_coded
        
    )
}
