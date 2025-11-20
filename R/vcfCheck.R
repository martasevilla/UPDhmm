#' Check variant quality, rename trio samples, and encode genotypes numerically
#' 
#' This function processes a VCF file and converts it into a *largeCollapsedVcf* 
#' object using the *VariantAnnotation* package. It renames the samples to standard 
#' names for subsequent analyses within the UPDhmm package and generates a numeric 
#' encoding of the trio genotypes for each variant, which is stored in the metadata 
#' column *geno_coded*.  
#' 
#' Additionally, the function includes an optional *quality_check* parameter. When 
#' enabled, it triggers warnings for variants that do not meet minimum quality 
#' thresholds based on read depth (RD) and genotype quality (GQ) values in the 
#' input VCF.
#'
#' @param largeCollapsedVcf A VCF object of class *CollapsedVCF*.
#' @param father Character string with the sample ID corresponding to the father.
#' @param mother Character string with the sample ID corresponding to the mother.
#' @param proband Character string with the sample ID corresponding to the proband.
#' @param check_quality Optional argument. TRUE/FALSE. If quality parameters 
#' want to be measured.
#' Default = FALSE
#'
#'
#' @return A *largeCollapsedVcf* object identical to the input with samples renamed 
#' to standard names for the trio and a new metadata column *geno_coded* containing 
#' the numeric encoding of the trio genotypes for each variant.
#' 
#' @examples
#' file <- system.file(package = "UPDhmm", "extdata", "test_het_mat.vcf.gz")
#' vcf  <- VariantAnnotation::readVcf(file)
#'
#' processedVcf <- vcfCheck(
#'     vcf,
#'     proband = "NA19675",
#'     mother  = "NA19678",
#'     father  = "NA19679"
#' )
#'
#' @export
#' 
vcfCheck <- function(
    largeCollapsedVcf,
    father,
    mother,
    proband,
    check_quality = FALSE) {
  
    # ---------------------------
    # Input validation
    # ---------------------------
    if (missing(largeCollapsedVcf)) stop("Argument 'largeCollapsedVcf' is missing.")
    if (!inherits(largeCollapsedVcf, "CollapsedVCF")) stop("Argument 'largeCollapsedVcf' must be a VCF object.")
    if (missing(proband)) stop("Argument 'proband' is missing.")
    if (missing(father)) stop("Argument 'father' is missing.")
    if (missing(mother)) stop("Argument 'mother' is missing.")
    if (!inherits(proband, "character")) stop("Argument 'proband' must be a character vector")
    if (!inherits(mother, "character")) stop("Argument 'mother' must be a character vector")
    if (!inherits(father, "character")) stop("Argument 'father' must be a character vector")
  
    # Extract genotype data
    geno_data <- VariantAnnotation::geno(largeCollapsedVcf)

    # Optional quality checks
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
    
    # Stop if any genotypes are not allowed
    if (!all(unique_genotypes %in% names(genotypes))) {
    invalid_genotypes <- unique_genotypes[!unique_genotypes %in% names(genotypes)]
    stop(paste("Error: The following genotypes are not valid:", 
               paste(unique(invalid_genotypes), collapse = ", ")))
    }

    # Save original sample IDs
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
