# create the expected output

file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
expected_vcf <- VariantAnnotation::readVcf(file)
colnames(expected_vcf) <- c("proband", "father", "mother")
S4Vectors::mcols(expected_vcf)$states <-
    c("iso_mat", "iso_mat", "iso_mat")


# read vcf
input <- VariantAnnotation::readVcf(file)
colnames(input) <- c("proband", "father", "mother")


utils::data("hmm")
hmm <- hmm
genotypes <- c(
    "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
    "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
)


test_that("Test if viterbi algorithm works", {
    out <- applyViterbi(
        largeCollapsedVcf = input,
        genotypes = genotypes,
        hmm = hmm
    )
    expect_s4_class(out, "CollapsedVCF")
    expect_equal(
        S4Vectors::mcols(out)$states,
        S4Vectors::mcols(expected_vcf)$states
    )
})
