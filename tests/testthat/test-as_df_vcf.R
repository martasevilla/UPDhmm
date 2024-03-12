expect_df <- data.frame(
  start = c(100017453, 100018844, 100144782),
  end = c(100017453, 100018844, 100144782),
  group = c("iso_fat", "iso_fat", "iso_fat"),
  seqnames = c("chr10", "chr10", "chr10"),
  genotype = c("313", "313", "313")
)

input <- VariantAnnotation::readVcf("test.vcf.gz")
colnames(input) <- c("proband", "father", "mother")
S4Vectors::mcols(input)$states <- c("iso_fat", "iso_fat", "iso_fat")

genotypes <-  c("0/0" = "1", "0/1" = "2","1/0" = "2", "1/1" = "3",
        "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3" )



test_that("Test if transofrmation of vcf with states to basic dataframe works",
          { out <- as_df_vcf(largecollapsedVcf = input,genotypes = genotypes)
  expect_s3_class(out, "data.frame")
  expect_equal(expect_df, out)
})
