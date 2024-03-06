expected_result <- data.frame(
  seqnames = "chr10",
  start = 100017453,
  end = 100144782,
  group = as.character("het_fat"),
  n_snps = 3,
  log_OR=31,
  p_value=2e-08
)



input <- VariantAnnotation::readVcf("test.vcf.gz")
colnames(input) <- c("proband", "father", "mother")
S4Vectors::mcols(input)$states <- c("het_fat", "het_fat", "het_fat")

position<-data.frame(start=100017453
             ,end=100144782
             ,group="het_fat"
             ,seqnames="chr10"
             ,n_snps=3)

utils::data("hmm")
hmm <- hmm
genotypes <-  c("0/0" = "1", "0/1" = "2","1/0" = "2", "1/1" = "3",
        "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3" )



test_that("Test if calculation of statistic parameters works", {
  out <- add_or(filtered_def_blocks_states = position,
        largecollapsedVcf = input,
        genotypes = genotypes,
        hmm=hmm)
  # round for errors in expect_equal
  out$log_OR <- floor(as.numeric(out$log_OR) * 10) / 10
  out$p_value <- floor(as.numeric(out$p_value) * 10^8) / 10^8
  out$seqnames<-as.character(out$seqnames)
  out$start<-as.numeric(out$start)
  out$end<-as.numeric(out$end)
  expect_equal(out, expected_result)
})
