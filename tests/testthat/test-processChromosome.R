# test-processChromosome.R

file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)

input <- vcfCheck(
  largeCollapsedVcf = input,
  father = "NA19689", mother = "NA19688",
  proband = "NA19685", check_quality = TRUE
)

# Split processed VCF by chromosome and select chromosome 6
split_vcf <- split(input, f = GenomicRanges::seqnames(input))
chr6 <- split_vcf[["6"]]

# Expected mean sequencing depth per individual for chromosome 6
total_mean_per_individual <- c(proband = 904/15, mother = 886/15, father = 902/15)

# Load the default HMM
utils::data("hmm", package = "UPDhmm", envir = environment())

# Identify which genotype codes correspond to Mendelian errors
emission_probs <- hmm$emissionProbs["normal", ]
mendelian_error_values <- names(emission_probs[emission_probs == min(emission_probs)])

# Expected output block 
expected_df <- data.frame(
  ID = "NA19685",
  seqnames = "6",
  start = 32489853,
  end=  33499925,
  group = "het_mat",
  n_snps = 5L,
  ratio_proband = 0.978982,
  ratio_mother = 1.002257,
  ratio_father = 0.951220,
  n_mendelian_error = 3L
)


# ------------------------------------------------------------------------- #
# Test processChromosome() with add_ratios = FALSE
# ------------------------------------------------------------------------- #

test_that("processChromosome works with valid chromosome input", {
  
  out <- processChromosome(vcf_chr = chr6, hmm = hmm, mendelian_error_values = mendelian_error_values)
  
  out$seqnames <- as.character(out$seqnames)  
  out <- as.data.frame(out)
  
  expected_no_ratio <- expected_df[, !(names(expected_df) %in% c("ratio_proband", "ratio_mother", "ratio_father"))]
  out_no_ratio <- out[, names(expected_no_ratio), drop = FALSE]
  
  expect_equal(out_no_ratio, expected_no_ratio)
  expect_s3_class(out, "data.frame")

})

# ------------------------------------------------------------------------- #
# Test processChromosome() with add_ratios = TRUE using DP field
# ------------------------------------------------------------------------- #
test_that("processChromosome works with valid chromosome input", {
  
  out <- processChromosome(vcf_chr = chr6, hmm = hmm, add_ratios = TRUE, field_DP = "DP", total_mean = total_mean_per_individual, mendelian_error_values = mendelian_error_values)
  
  out$seqnames <- as.character(out$seqnames)
  out <- as.data.frame(out)
  
  expect_equal(out, expected_df, tolerance = 1e-6)
  expect_s3_class(out, "data.frame")
  
})

# ------------------------------------------------------------------------- #
# Test processChromosome() with add_ratios = TRUE using AD field
# ------------------------------------------------------------------------- #
test_that("processChromosome works with valid chromosome input", {
  
  out <- processChromosome(vcf_chr = chr6, hmm = hmm, add_ratios = TRUE, field_DP = "AD", total_mean = total_mean_per_individual, mendelian_error_values = mendelian_error_values)
  
  out$seqnames <- as.character(out$seqnames)
  out <- as.data.frame(out)
  
  expect_equal(out, expected_df, tolerance = 1e-6)
  expect_s3_class(out, "data.frame")
  
})

# ------------------------------------------------------------------------- #
# Test processChromosome() with add_ratios = TRUE and no field_DP defined
# ------------------------------------------------------------------------- #
test_that("processChromosome works with valid chromosome input", {
  
  out <- processChromosome(vcf_chr = chr6, hmm = hmm, add_ratios = TRUE, total_mean = total_mean_per_individual, mendelian_error_values = mendelian_error_values)
  out$seqnames <- as.character(out$seqnames)
  out <- as.data.frame(out)
  
  expect_equal(out, expected_df, tolerance = 1e-6)
  expect_s3_class(out, "data.frame")
  
})

# ------------------------------------------------------------------------- #
# Modify DP and AD to introduce NA values
# ------------------------------------------------------------------------- #
g_dp <- VariantAnnotation::geno(chr6)$DP
g_ad <- VariantAnnotation::geno(chr6)$AD

# Introduce NA in proband DP at the first variant
g_dp[1, "proband"] <- NA
VariantAnnotation::geno(chr6)$DP <- g_dp

# Introduce NA in proband AD (both alleles NA) at the first variant
g_ad[1, "proband"][[1]] <- c(NA,NA)
VariantAnnotation::geno(chr6)$AD <- g_ad

# Expected proband ratio after introducing NA values
expected_df$ratio_proband <- 0.974526

total_mean_per_individual <- c(proband = 844/14, mother = 886/15, father = 902/15)

# ------------------------------------------------------------------------- #
# Repeat the three tests under conditions where NA values are present
# ------------------------------------------------------------------------- #

test_that("processChromosome works with valid chromosome input", {
  
  out <- processChromosome(vcf_chr = chr6, hmm = hmm, add_ratios = TRUE, field_DP = "DP", total_mean = total_mean_per_individual, mendelian_error_values = mendelian_error_values)
  
  out$seqnames <- as.character(out$seqnames)
  out <- as.data.frame(out)
  
  expect_equal(out, expected_df, tolerance = 1e-6)
  expect_s3_class(out, "data.frame")
  
})


test_that("processChromosome works with valid chromosome input", {
  
  out <- processChromosome(vcf_chr = chr6, hmm = hmm, add_ratios = TRUE, field_DP = "AD", total_mean = total_mean_per_individual, mendelian_error_values = mendelian_error_values)
  
  out$seqnames <- as.character(out$seqnames)
  out <- as.data.frame(out)
  
  expect_equal(out, expected_df, tolerance = 1e-6)
  expect_s3_class(out, "data.frame")
  
})

test_that("processChromosome works with valid chromosome input", {
  
  out <- processChromosome(vcf_chr = chr6, hmm = hmm, add_ratios = TRUE, total_mean = total_mean_per_individual, mendelian_error_values = mendelian_error_values)
  
  out$seqnames <- as.character(out$seqnames)
  out <- as.data.frame(out)
  
  expect_equal(out, expected_df, tolerance = 1e-6)
  expect_s3_class(out, "data.frame")
  
})
