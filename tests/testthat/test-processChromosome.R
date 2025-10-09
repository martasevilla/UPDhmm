# test-processChromosome.R
file <- system.file(package = "UPDhmm", "extdata", "test_het_mat.vcf.gz")
vcf <- VariantAnnotation::readVcf(file)
processedVcf <- vcfCheck(
  vcf,
  proband = "NA19675", mother = "NA19678", father = "NA19679"
)
# Split by chromosome
split_vcf <- split(processedVcf, f = GenomicRanges::seqnames(vcf))
chr1 <- split_vcf[[1]]  # pick first chromosome subset

# Load default HMM
utils::data("hmm", package = "UPDhmm", envir = environment())

genotypes <- c(
  "0/0" = "1", "0/1" = "2", "1/0" = "2", "1/1" = "3",
  "0|0" = "1", "0|1" = "2", "1|0" = "2", "1|1" = "3"
)

expected_df <- data.frame(
  ID = "NA19675",
  start = 1254841,
  end= 248845499,
  group = "normal",
  seqnames = "1",
  n_snps = 340
)



test_that("processChromosome works with valid chromosome input", {
  
  out <- processChromosome(chr1, hmm = hmm, genotypes = genotypes)
  expect_equal(out, expected_df)

})

