expected_result <- data.frame(
  seqnames = "6",
  start = 32489853,
  end = 32489876,
  group = "iso_mat",
  n_snps = 3,
  log_likelihood = 21.1,
  p_value = 4e-6,
  n_mendelian_error=2
)


file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)
colnames(input) <- c("proband", "father", "mother")
S4Vectors::mcols(input)$states <- c("iso_mat", "iso_mat", "iso_mat")

position <- data.frame(
  start = 32489853,
  end = 32489876,
  group = "iso_mat",
  seqnames = "6",
  n_snps = 3
)

utils::data("hmm")
hmm <- hmm
genotypes <- c(
  "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
  "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
)



test_that("Test if calculation of statistic parameters works", {
  out <- addOr(
    filtered_def_blocks_states = position,
    largeCollapsedVcf = input,
    genotypes = genotypes,
    hmm = hmm
  )
  # round for errors in expect_equal
  out$log_likelihood <- floor(as.numeric(out$log_likelihood) * 10) / 10
  out$p_value <- floor(as.numeric(out$p_value) * 10^6) / 10^6
  out$seqnames <- as.character(out$seqnames)
  out$start <- as.numeric(out$start)
  out$end <- as.numeric(out$end)
  expect_equal(out, expected_result)
})
