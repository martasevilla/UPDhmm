#' Check quality parameters (optional) and change IDs.
#'
#' This function takes a VCF file and converts it into a largeCollapsedVcf
#' object using the VariantAnnotation package. It also rename the sample for 
#' subsequent steps needed in UPDhmm package.
#' Additionally, it features an optional parameter, quality_check, which triggers warnings 
#' when variants lack sufficient quality based on RD and GQ parameters in the input VCF.
#'
#' @param largeCollapsedVcf The file in largeCollapsedVcf format.
#' @param father Name of the father's sample.
#' @param mother Name of the mother's sample.
#' @param proband Name of the proband's sample.
#' @param check_quality Optional argument. TRUE/FALSE. If quality parameters 
#' want to be measured.
#' Default = FALSE
#'
#' @return largeCollapsedVcf (VariantAnnotation VCF format).
#' @export
#' @examples
#' fl <- system.file("extdata", "test_het_mat.vcf.gz", package = "UPDhmm")
#' vcf <- VariantAnnotation::readVcf(fl)
#' processedVcf <-
#'    vcfCheck(vcf, proband = "Sample1", mother = "Sample3", father = "Sample2")
#'
vcfCheck <- function(
    largeCollapsedVcf,
    father,
    mother,
    proband,
    check_quality = FALSE) {
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
    # Extract genotype data
    geno_data <- VariantAnnotation::geno(largeCollapsedVcf)
    
    # Quality parameters
    if (isTRUE(check_quality)) {
      if (any(geno_data$GQ < 20 | is.na(geno_data$GQ)) == TRUE) {
        message("No filter quality (GQ) parameter used")
      }
      if (any(geno_data$DP < 30 | is.na(geno_data$DP)) == TRUE) {
        message("No filter quality (RD) parameter used")
      }
    }
    
    # Define allowed genotypes and numeric codes
    genotypes <- c(
      "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
      "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
    )
    
    unique_genotypes <- unique(geno_data$GT)
    
    if (!all(unique_genotypes %in% names(genotypes))) {
      invalid_genotypes <- unique_genotypes[!unique_genotypes %in% names(genotypes)]
      stop(paste("Error: The following genotypes are not valid:", 
                 paste(unique(invalid_genotypes), collapse = ", ")))
    }
    
    # Save original IDs in colData
    SummarizedExperiment::colData(largeCollapsedVcf)$ID <- colnames(largeCollapsedVcf)
    
    # Change names in vcf for subsequent steps
    colnames(largeCollapsedVcf)[colnames(largeCollapsedVcf) 
                                == father] <- "father"
    colnames(largeCollapsedVcf)[colnames(largeCollapsedVcf) 
                                == mother] <- "mother"
    colnames(largeCollapsedVcf)[colnames(largeCollapsedVcf) 
                                == proband] <- "proband"
    
    # Generate numeric encoding of trio genotypes
    geno_uncoded <- VariantAnnotation::geno(largeCollapsedVcf)$GT
    geno_coded <- paste0(
      genotypes[geno_uncoded[, "father"]],
      genotypes[geno_uncoded[, "mother"]],
      genotypes[geno_uncoded[, "proband"]]
    )
    # Each string is 3 characters: father + mother + proband; values: 1=hom_ref, 2=het, 3=hom_alt
    
    # Store the numeric genotype string in the VCF metadata column 'geno_coded'
    S4Vectors::mcols(largeCollapsedVcf)$geno_coded <- geno_coded

    
    return(largeCollapsedVcf)
    
}