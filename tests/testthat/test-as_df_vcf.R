expect_df <- data.frame(
  start = c(100017453, 100018844, 100144782),
  end = c(100017453, 100018844, 100144782),
  group = c("het_fat", "het_fat", "het_fat"),
  seqnames = c("chr10", "chr10", "chr10"),
  genotype = c("313", "313", "313")
)

input <- VariantAnnotation::readVcf("test.vcf.gz")
colnames(input) <- c("proband", "father", "mother")
GenomicRanges::elementMetadata(input)$metadata <- c("het_fat", "het_fat", "het_fat")


test_that("Test if transofrmation of vcf with states to basic dataframe works", {
  out <- as_df_vcf(input)

  expect_s3_class(out, "data.frame")
  expect_equal(expect_df, out)
})
