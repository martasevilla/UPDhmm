expected_df <- data.frame(
  start = 100017453,
  end = 100144782,
  group = as.character("het_fat"),
  n_snps = 3,
  seqnames = "chr10",
  row.names = "1"
)

expected_df$group <- as.character(expected_df$group)


input <- data.frame(
  start = c(100017453, 100018844, 100144782),
  end = c(100017453, 100018844, 100144782),
  group = c("het_fat", "het_fat", "het_fat"),
  seqnames = c("chr10", "chr10", "chr10")
)
test_that("Check if simplification into blocks is working", {
  out <- blocks_vcf(input)

  expect_equal(out, expected_df)
  expect_s3_class(out, "data.frame")
})
