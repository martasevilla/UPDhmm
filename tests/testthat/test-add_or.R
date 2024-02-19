expected_result <- c(log_OR = 31, p_value = 2e-08)


split_vcf_df <- list(
  chr10 =
    data.frame(
      start = c(100017453, 100018844, 100144782),
      end = c(100017453, 100018844, 100144782),
      group = c("het_fat", "het_fat", "het_fat"),
      seqnames = c("chr10", "chr10", "chr10"),
      genotype = c("313", "313", "313")
    )
)


test_that("Test if calculation of statistic parameters works", {
  out <- add_or(100017453, 100144782, "chr10", "het_fat", split_vcf_df)
  # round for errors in expect_equal
  out["log_OR"] <- floor(as.numeric(out["log_OR"]) * 10) / 10
  out["p_value"] <- floor(as.numeric(out["p_value"]) * 10^8) / 10^8

  expect_equal(out, expected_result)
})
