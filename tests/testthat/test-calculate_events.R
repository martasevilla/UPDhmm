expected_def_blocks <- data.frame(
  seqnames = "chr10",
  start = 100017453,
  end = 100144782,
  group = "het_fat",
  n_snps = 3,
  log_OR = 31,
  p_value = 2e-08
)


input <- VariantAnnotation::readVcf("test.vcf.gz")
colnames(input) <- c("proband", "father", "mother")
test_that("Test if the general function works", {
  out <- calculate_events(input)
  out$seqnames<-as.character(out$seqnames)
  out$log_OR <- floor(as.numeric(out$log_OR) * 10) / 10
  out$p_value <- floor(as.numeric(out$p_value) * 10^8) / 10^8
  out<-as.data.frame(out)
  expect_equal(expected_def_blocks, out)
  expect_s3_class(out, "data.frame")
})
