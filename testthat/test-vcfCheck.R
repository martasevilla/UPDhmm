file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
expected_vcf <- VariantAnnotation::readVcf(file)
colnames(expected_vcf) <- c("proband", "father", "mother")

test_that("Test if the vcf loading works", {
    input_vcf <- VariantAnnotation::readVcf(file)
    input <- vcfCheck(input_vcf,
        father = "NA19679", mother = "NA19678",
        proband = "NA19675", check_quality = TRUE
    )
    expect_s4_class(input, "CollapsedVCF")
    expect_equal(input, expected_vcf)
})
