#' check quality parameters and change ids
#'
#' This function will take the vcf (must be gz with an index file and will
#' convert it into a large collapsedVcf object using variantanntoation package)
#' it will have one additional parameter, quality_check, that will give warnigns
#' if there are positions with no good quality.
#'
#' @param path_vcf The file in largeCollapsed format
#' @param check_quality If desired, quality parameters can be measured
#' @param father Introduce father's name sample
#' @param mother Introduce mother's name sample
#' @param proband Introduce proband's name sample
#'
#' @return LargecollapdsedVCf (VariantAnnotation vcf format)
#' @export
#' @examples
#' fl <- system.file("extdata", "test.vcf.gz", package = "UPDhmm")
#' vcf <- VariantAnnotation::readVcf(fl)
#' vcf_check(vcf, proband = "Sample1", mother = "Sample3", father = "Sample2")
#'
#'
#'
#'
#'
#'
# sacar esto de fuera
vcf_check <- function(path_vcf, check_quality = NULL, father = NULL,
                      mother = NULL, proband = NULL) {
  if (isTRUE(check_quality)) {
    # quality parameters
    if (any(VariantAnnotation::geno(path_vcf)$GQ < 20 |
      is.na(VariantAnnotation::geno(path_vcf)$GQ)) == TRUE) {
      warning("No filter quality (GQ) parameter used")
    }
    if (any(VariantAnnotation::geno(path_vcf)$DP < 30 |
      is.na(VariantAnnotation::geno(path_vcf)$DP)) == TRUE) {
      warning("No filter quality (RD) parameter used")
    }


    # Update sample names in the header lines
    colnames(path_vcf) <- gsub(father, "father", colnames(path_vcf))
    colnames(path_vcf) <- gsub(mother, "mother", colnames(path_vcf))
    colnames(path_vcf) <- gsub(proband, "proband", colnames(path_vcf))
  } else {
    # Update sample names in the header lines
    colnames(path_vcf) <- gsub(father, "father", colnames(path_vcf))
    colnames(path_vcf) <- gsub(mother, "mother", colnames(path_vcf))
    colnames(path_vcf) <- gsub(proband, "proband", colnames(path_vcf))
  }

  return(path_vcf)
}
