#' Check quality parameters (optional) and change IDs.
#'
#' This function takes a VCF file and converts it into a large collapsedVcf
#' object using the VariantAnnotation package. It includes an optional parameter,
#' quality_check, which issues warnings if positions lack good quality based on
#' RD and GQ parameters in the input VCF.
#'
#' @param path_vcf The file in largeCollapsed format.
#' @param check_quality TRUE/FALSE. If desired, quality parameters can be measured.
#' @param father Name of the father's sample.
#' @param mother Name of the mother's sample.
#' @param proband Name of the proband's sample.
#'
#' @return LargecollapsedVCF (VariantAnnotation VCF format).
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
vcf_check <- function(path_vcf, check_quality = NULL, father = NULL,
                      mother = NULL, proband = NULL) {
  # Quality parameters
  if (isTRUE(check_quality)) {
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
  }
  else {
    # No quality parameters, just update sample names in the header lines
    colnames(path_vcf) <- gsub(father, "father", colnames(path_vcf))
    colnames(path_vcf) <- gsub(mother, "mother", colnames(path_vcf))
    colnames(path_vcf) <- gsub(proband, "proband", colnames(path_vcf))
  }

  return(path_vcf)
}
