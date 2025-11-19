# test-vcfCheck.R
file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
expected_vcf <- VariantAnnotation::readVcf(file)
SummarizedExperiment::colData(expected_vcf)$ID <- colnames(expected_vcf)
colnames(expected_vcf) <- c("proband", "father", "mother")

# Expected numeric genotype encoding for chromosome 6
expected_geno_coded <- c("133", "133", "121", "122", "133")

test_that("Test if the vcf loading works", {
  input_vcf <- VariantAnnotation::readVcf(file)
  input <- vcfCheck(input_vcf,
                    father = "NA19689", mother = "NA19688",
                    proband = "NA19685", check_quality = TRUE
  )
  
  # Check that the returned object is of the correct S4 class
  expect_s4_class(input, "CollapsedVCF")
  
  # Verify that sample names were renamed correctly
  expect_equal(colnames(input), c("proband", "father", "mother"))
  
  # Confirm that the new metadata column 'geno_coded' was created
  expect_true("geno_coded" %in% colnames(S4Vectors::mcols(input)))
  
  # Split the VCF by chromosome and extract chromosome 6
  split_input <- split(input, f = GenomicRanges::seqnames(input))
  chr6 <- split_input[["6"]]
  
  # Retrieve the computed genotype encoding
  geno_coded <- S4Vectors::mcols(chr6)$geno_coded
  
  # Compare with expected output
  expect_equal(geno_coded, expected_geno_coded)
})

# ------------------------------------------------------------------------- #

test_that("vcfCheck reports low read depth", {
  input_vcf <- VariantAnnotation::readVcf(file)
  geno_data <- VariantAnnotation::geno(input_vcf)
  
  # Manually lower DP for the first variant to trigger the warning message
  geno_data$DP[1, ] <- 10
  VariantAnnotation::geno(input_vcf) <- geno_data
  
  # Expect a message about insufficient read depth
  input <- expect_message(
    vcfCheck(input_vcf,
             father = "NA19689", mother = "NA19688",
             proband = "NA19685", check_quality = TRUE),
    "No filter quality \\(RD\\) parameter used"
  )
  
  # Check that the returned object is still valid
  expect_s4_class(input, "CollapsedVCF")
  
  # Ensure the sample columns were renamed as expected
  expect_equal(colnames(input), c("proband", "father", "mother"))
  
  # Verify that geno_coded was generated
  expect_true("geno_coded" %in% colnames(S4Vectors::mcols(input)))
})
