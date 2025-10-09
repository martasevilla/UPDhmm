file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
expected_vcf <- VariantAnnotation::readVcf(file)
SummarizedExperiment::colData(expected_vcf)$ID <- colnames(expected_vcf)
colnames(expected_vcf) <- c("proband", "father", "mother")

test_that("Test if the vcf loading works", {
    input_vcf <- VariantAnnotation::readVcf(file)
    input <- vcfCheck(input_vcf,
        father = "NA19689", mother = "NA19688",
        proband = "NA19685", check_quality = TRUE
    )
    expect_s4_class(input, "CollapsedVCF")
    expect_equal(input, expected_vcf)
})
