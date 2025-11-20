# test-applyViterbi.R

file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
expected_vcf <- VariantAnnotation::readVcf(file)
colnames(expected_vcf) <- c("proband", "father", "mother")
S4Vectors::mcols(expected_vcf)$states <-
    c("iso_mat", "iso_mat", "iso_mat")


input <- VariantAnnotation::readVcf(file)

input <- vcfCheck(
  largeCollapsedVcf = input,
  father = "NA19689", mother = "NA19688",
  proband = "NA19685", check_quality = TRUE
)

# Load the default HMM
utils::data("hmm", package = "UPDhmm")

test_that("Test if viterbi algorithm works", {
    out <- applyViterbi(
        largeCollapsedVcf = input,
        hmm = hmm
    )
    expect_s4_class(out, "CollapsedVCF")
    expect_equal(
        S4Vectors::mcols(out)$states,
        S4Vectors::mcols(expected_vcf)$states
    )
    expect_equal(length(S4Vectors::mcols(out)$states), nrow(out))
    
})
