expect_df <- data.frame(
    start = c(32489853, 32489856, 32489876),
    end = c(32489853, 32489856, 32489876),
    group = c("iso_mat", "iso_mat", "iso_mat"),
    seqnames = c("6", "6", "6"),
    genotype = c("133", "133", "121")
)


file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)
colnames(input) <- c("proband", "father", "mother")
S4Vectors::mcols(input)$states <- c("iso_mat", "iso_mat", "iso_mat")

genotypes <- c(
    "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
    "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
)



test_that("Test if transformation of vcf with states to basic dataframe works", {
    out <- asDfVcf(largeCollapsedVcf = input, genotypes = genotypes)
    expect_s3_class(out, "data.frame")
    expect_equal(expect_df, out)
})
