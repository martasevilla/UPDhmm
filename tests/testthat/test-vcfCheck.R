expected_vcf <- VariantAnnotation::readVcf("test.vcf.gz")
colnames(expected_vcf) <- c("proband", "father", "mother")

test_that("Test if the vcf loading works", {
    input_vcf <- VariantAnnotation::readVcf("test.vcf.gz")
    input <- vcfCheck(input_vcf,
    father = "Sample2", mother = "Sample3",
    proband = "Sample1", check_quality = TRUE)
    expect_s4_class(input, "CollapsedVCF")
    expect_equal(input, expected_vcf)
})
