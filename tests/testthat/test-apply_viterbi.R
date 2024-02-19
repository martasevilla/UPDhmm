# create the expected output
expected_vcf <- VariantAnnotation::readVcf("test.vcf.gz")
colnames(expected_vcf) <- c("proband", "father", "mother")
GenomicRanges::elementMetadata(expected_vcf)$metadata <-
  c("het_fat", "het_fat", "het_fat")



# read vcf
input <- VariantAnnotation::readVcf("test.vcf.gz")
colnames(input) <- c("proband", "father", "mother")

test_that("Check if viterbi algorithm is working", {
  #### no meter otras funciones del paquete

  out <- apply_viterbi(input)
  expect_s4_class(out, "CollapsedVCF")
  expect_equal(
    GenomicRanges::elementMetadata(out)$metadata,
    GenomicRanges::elementMetadata(expected_vcf)$metadata
  )
})
