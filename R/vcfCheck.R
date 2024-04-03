#' Check quality parameters (optional) and change IDs.
#'
<<<<<<< HEAD
#' This function takes a VCF file and converts it into a largeCollapsedVcf
=======
#' This function takes a VCF file and converts it into a large collapsedVcf
>>>>>>> upstream/devel
#' object using the VariantAnnotation package. It includes an optional
#' parameter, quality_check, which issues warnings if positions lack good
#' quality based on RD and GQ parameters in the input VCF.
#'
<<<<<<< HEAD
#' @param largeCollapsedVcf The file in largeCollapsedVcf format.
=======
#' @param largecollapsedVcf The file in largeCollapsed format.
>>>>>>> upstream/devel
#' @param check_quality TRUE/FALSE. If quality parameters want to be measured.
#' Default = FALSE
#' @param father Name of the father's sample.
#' @param mother Name of the mother's sample.
#' @param proband Name of the proband's sample.
#'
<<<<<<< HEAD
#' @return largeCollapsedVcf (VariantAnnotation VCF format).
=======
#' @return LargecollapsedVCF (VariantAnnotation VCF format).
>>>>>>> upstream/devel
#' @export
#' @examples
#' fl <- system.file("extdata", "test_het_mat.vcf.gz", package = "UPDhmm")
#' vcf <- VariantAnnotation::readVcf(fl)
<<<<<<< HEAD
#' processedVcf <-
#'    vcfCheck(vcf, proband = "Sample1", mother = "Sample3", father = "Sample2")
#'
vcfCheck <- function(
    largeCollapsedVcf,
    check_quality = FALSE,
    father,
    mother,
    proband) {
    # Check if `largeCollapsedVcf` is provided

    if (missing(largeCollapsedVcf)) {
        stop("Argument 'largeCollapsedVcf' is missing.")
    }

    # Check if `largeCollapsedVcf` is a VCF object
    if (!inherits(largeCollapsedVcf, "CollapsedVCF")) {
        stop("Argument 'largeCollapsedVcf' must be a VCF object.")
    }

    # Check if `proband`,`father` and `mother` is provided

    if (missing(proband)) {
        stop("Argument 'proband' is missing.")
    }
    if (missing(father)) {
        stop("Argument 'father' is missing.")
    }

    if (missing(mother)) {
        stop("Argument 'mother' is missing.")
    }


    # Check if `largeCollapsedVcf` is a VCF object
    if (!inherits(largeCollapsedVcf, "CollapsedVCF")) {
        stop("Argument 'largeCollapsedVcf' must be a VCF object.")
    }

    # Check if `proband` ,`father` and `mother` is a character vector
    if (!inherits(proband, "character")) {
        stop("Argument 'proband' must be a character vector")
    }

    if (!inherits(mother, "character")) {
        stop("Argument 'mother' must be a character vector")
    }

    if (!inherits(father, "character")) {
        stop("Argument 'father' must be a character vector")
    }

    # Quality parameters
    if (isTRUE(check_quality)) {
        if (any(VariantAnnotation::geno(largeCollapsedVcf)$GQ < 20 |
            is.na(VariantAnnotation::geno(largeCollapsedVcf)$GQ)) == TRUE) {
            message("No filter quality (GQ) parameter used")
        }
        if (any(VariantAnnotation::geno(largeCollapsedVcf)$DP < 30 |
            is.na(VariantAnnotation::geno(largeCollapsedVcf)$DP)) == TRUE) {
            message("No filter quality (RD) parameter used")
        }
    }

    # Change names in vcf for subsequent steps

  colnames(largeCollapsedVcf)[colnames(largeCollapsedVcf) 
                              == father] <- "father"
  colnames(largeCollapsedVcf)[colnames(largeCollapsedVcf) 
                              == mother] <- "mother"
  colnames(largeCollapsedVcf)[colnames(largeCollapsedVcf) 
                              == proband] <- "proband"

    return(largeCollapsedVcf)
}
=======
#' processedVcf<-
#' vcfCheck(vcf, proband = "Sample1", mother = "Sample3", father = "Sample2")
#'
#'
#'
#'
#'
#'
vcfCheck <- function(
largecollapsedVcf,
check_quality = FALSE,
father = NULL,
mother = NULL,
proband = NULL) {
  # Quality parameters
if (isTRUE(check_quality)) {
if (any(VariantAnnotation::geno(largecollapsedVcf)$GQ < 20 |
is.na(VariantAnnotation::geno(largecollapsedVcf)$GQ)) == TRUE) {
message("No filter quality (GQ) parameter used")
}
if (any(VariantAnnotation::geno(largecollapsedVcf)$DP < 30 |
is.na(VariantAnnotation::geno(largecollapsedVcf)$DP)) == TRUE) {
message("No filter quality (RD) parameter used")
}}

#Change names in vcf for subsequent steps
colnames(largecollapsedVcf) <-
gsub(father, "father",colnames(largecollapsedVcf))
colnames(largecollapsedVcf) <-
gsub(mother, "mother",colnames(largecollapsedVcf))
colnames(largecollapsedVcf) <-
gsub(proband, "proband",colnames(largecollapsedVcf))
return(largecollapsedVcf)

}

>>>>>>> upstream/devel
