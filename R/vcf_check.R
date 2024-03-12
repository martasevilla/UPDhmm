#' Check quality parameters (optional) and change IDs.
#'
#' This function takes a VCF file and converts it into a large collapsedVcf
#' object using the VariantAnnotation package. It includes an optional
#' parameter, quality_check, which issues warnings if positions lack good
#' quality based on RD and GQ parameters in the input VCF.
#'
#' @param largecollapsedVcf The file in largeCollapsed format.
#' @param check_quality TRUE/FALSE. If quality parameters want to be measured.
#' Default = FALSE
#' @param father Name of the father's sample.
#' @param mother Name of the mother's sample.
#' @param proband Name of the proband's sample.
#'
#' @return LargecollapsedVCF (VariantAnnotation VCF format).
#' @export
#' @examples
#' fl <- system.file("extdata", "test_het_mat.vcf.gz", package = "UPDhmm")
#' vcf <- VariantAnnotation::readVcf(fl)
#' vcf_check(vcf, proband = "Sample1", mother = "Sample3", father = "Sample2")
#'
#'
#'
#'
#'
#'
vcf_check <- function(largecollapsedVcf, check_quality = FALSE, father = NULL,
            mother = NULL, proband = NULL) {
  # Quality parameters
  if (isTRUE(check_quality)) {
  if (any(VariantAnnotation::geno(largecollapsedVcf)$GQ < 20 |
    is.na(VariantAnnotation::geno(largecollapsedVcf)$GQ)) == TRUE) {
    warning("No filter quality (GQ) parameter used")
  }
  if (any(VariantAnnotation::geno(largecollapsedVcf)$DP < 30 |
    is.na(VariantAnnotation::geno(largecollapsedVcf)$DP)) == TRUE) {
    warning("No filter quality (RD) parameter used")
  }

  }

#Change names in vcf for subsequent steps
  colnames(largecollapsedVcf) <- gsub(father, "father",
                    colnames(largecollapsedVcf))
  colnames(largecollapsedVcf) <- gsub(mother, "mother",
                    colnames(largecollapsedVcf))
  colnames(largecollapsedVcf) <- gsub(proband, "proband",
                    colnames(largecollapsedVcf))
  return(largecollapsedVcf)
}

