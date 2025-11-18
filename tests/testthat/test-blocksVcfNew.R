# test-blocksVcfNew.R

# Expected output block 
expected_df <- data.frame(
  ID = "NA19685",
  seqnames = "6",
  start = 32489853,
  end=  33499925,
  group = "iso_mat",
  n_snps = 5L,
  geno_coded = I(setNames(list(c("133", "133", "121", "122", "133")), "1")),
  ratio_proband = 0.978982,
  ratio_mother = 1.002257,
  ratio_father = 0.951220
)

file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)

input <- vcfCheck(
  largeCollapsedVcf = input,
  father = "NA19689", mother = "NA19688",
  proband = "NA19685", check_quality = TRUE
)

# Load the default HMM
utils::data("hmm", package = "UPDhmm")
hmm <- hmm
input <- applyViterbi(input, hmm)

# Split processed VCF by chromosome and select chromosome 6
split_vcf <- split(input, f = GenomicRanges::seqnames(input))
chr6 <- split_vcf[[6]]

# Expected mean sequencing depth per individual for chromosome 6
total_mean <- c(proband = 904/15, mother = 886/15, father = 902/15)

# ------------------------------------------------------------------------- #
# Test blocksVcfNew() with add_ratios = FALSE
# ------------------------------------------------------------------------- #

test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(largeCollapsedVcf = chr6)

    out <- as.data.frame(out)
    
    expected_no_ratio <- expected_df[, !(names(expected_df) %in% c("ratio_proband", "ratio_mother", "ratio_father"))]
    
    expect_equal(out, expected_no_ratio)
    expect_s3_class(out, "data.frame")
})

# ------------------------------------------------------------------------- #
# Test blocksVcfNew() with add_ratios = TRUE using DP field
# ------------------------------------------------------------------------- #

test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(largeCollapsedVcf = chr6, add_ratios = TRUE, field_DP = "DP", total_mean = total_mean)

    out <- as.data.frame(out)
    
    out$ratio_proband <- round(out$ratio_proband, 6)
    out$ratio_mother <- round(out$ratio_mother, 6)
    out$ratio_father <- round(out$ratio_father, 6)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})

# ------------------------------------------------------------------------- #
# Test blocksVcfNew() with add_ratios = TRUE using AD field
# ------------------------------------------------------------------------- #

test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(largeCollapsedVcf = chr6, add_ratios = TRUE, field_DP = "AD", total_mean = total_mean)

    out <- as.data.frame(out)
    
    out$ratio_proband <- round(out$ratio_proband, 6)
    out$ratio_mother <- round(out$ratio_mother, 6)
    out$ratio_father <- round(out$ratio_father, 6)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})

# ------------------------------------------------------------------------- #
# Test blocksVcfNew() with add_ratios = TRUE and no field_DP defined
# ------------------------------------------------------------------------- #

test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(largeCollapsedVcf = chr6, add_ratios = TRUE, total_mean = total_mean)

    out <- as.data.frame(out)
    
    out$ratio_proband <- round(out$ratio_proband, 6)
    out$ratio_mother <- round(out$ratio_mother, 6)
    out$ratio_father <- round(out$ratio_father, 6)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})

# ------------------------------------------------------------------------- #
# Modify DP and AD to introduce NA values
# ------------------------------------------------------------------------- #
g_dp <- VariantAnnotation::geno(input)$DP
g_ad <- VariantAnnotation::geno(input)$AD

# Introduce NA in proband DP at the first variant
g_dp[1, "proband"] <- NA
VariantAnnotation::geno(input)$DP <- g_dp

# Introduce NA in proband AD (both alleles NA) at the first variant
g_ad[1, "proband"][[1]] <- c(NA,NA)
VariantAnnotation::geno(input)$AD <- g_ad

split_vcf <- split(input, f = GenomicRanges::seqnames(input))
chr6 <- split_vcf[[6]]

# Expected proband ratio after introducing NA values
expected_df$ratio_proband <- 0.974526

total_mean <- c(proband = 844/14, mother = 886/15, father = 902/15)


# ------------------------------------------------------------------------- #
# Repeat the three tests under conditions where NA values are present
# ------------------------------------------------------------------------- #

test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(largeCollapsedVcf = chr6, add_ratios = TRUE, field_DP = "DP", total_mean = total_mean)

    out <- as.data.frame(out)
    
    out$ratio_proband <- round(out$ratio_proband, 6)
    out$ratio_mother <- round(out$ratio_mother, 6)
    out$ratio_father <- round(out$ratio_father, 6)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})


test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(largeCollapsedVcf = chr6, add_ratios = TRUE, field_DP = "AD", total_mean = total_mean)

    out <- as.data.frame(out)
    
    out$ratio_proband <- round(out$ratio_proband, 6)
    out$ratio_mother <- round(out$ratio_mother, 6)
    out$ratio_father <- round(out$ratio_father, 6)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})

test_that("Test if simplification into blocks works", {
    out <- blocksVcfNew(largeCollapsedVcf = chr6, add_ratios = TRUE, total_mean = total_mean)

    out <- as.data.frame(out)
    
    out$ratio_proband <- round(out$ratio_proband, 6)
    out$ratio_mother <- round(out$ratio_mother, 6)
    out$ratio_father <- round(out$ratio_father, 6)

    expect_equal(out, expected_df)
    expect_s3_class(out, "data.frame")
})
