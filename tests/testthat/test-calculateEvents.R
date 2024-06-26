expected_def_blocks <- data.frame(
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
test_that("Test if the general function works", {
    out <- calculateEvents(input)
    out$seqnames <- as.character(out$seqnames)
    out$log_likelihood <- floor(as.numeric(out$log_likelihood) * 10) / 10
    out$p_value <- floor(as.numeric(out$p_value) * 10^6) / 10^6
    out <- as.data.frame(out)
    expect_equal(expected_def_blocks, out)
    expect_s3_class(out, "data.frame")
})





new_hmm<-list(
  States = c("normal", "iso_fat", "iso_mat", "het_fat", "het_mat")
  ,
  Symbols = c("111", "112", "113", "121", "122", "123", "131", "132", "133", 
              "211", "212", "213", "221", "222", "223", "231", "232", "233", 
              "311", "312", "313", "321", "322", "323", "331", "332", "333", 
              "000")
  ,
  startProbs = c(normal = 0.996, iso_fat = 0.001, iso_mat = 0.001, 
                 het_fat = 0.001, het_mat = 0.001)
  ,
  transProbs = matrix(
    c(
      0.99996, 0.00001, 0.00001, 0.00001, 0.00001,
      0.00001, 0.99996, 0.00001, 0.00001, 0.00001,
      0.00001, 0.00001, 0.99996, 0.00001, 0.00001,
      0.00001, 0.00001, 0.00001, 0.99996, 0.00001,
      0.00001, 0.00001, 0.00001, 0.00001, 0.99996
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      from = c("normal", "iso_fat", "iso_mat", "het_fat", "het_mat"),
      to = c("normal", "iso_fat", "iso_mat", "het_fat", "het_mat")
    )
  )
  
  
  ,
  emissionProbs = matrix(
    c(0.09250,0.00001,0.00001,0.09250,0.09250,0.00001,0.00001,0.09250,0.00001,0.09250,0.09250,0.00001,0.09250,0.12500,0.09250,0.00001,0.09250,0.09250,0.00001,0.09250,0.00001,0.00001,0.09250,0.09250,0.00001,0.00001,0.09250,0.00006,0.09250,0.00001,0.00001,0.12500,0.00001,0.00001,0.09250,0.00001,0.00001,0.09250,0.00001,0.09250,0.12500,0.00001,0.12500,0.09250,0.00001,0.09250,0.00001,0.00001,0.09250,0.00001,0.00001,0.12500,0.00001,0.00001,0.09250,0.00003,0.09250,0.00001,0.00001,0.09250,0.00001,0.09250,0.00001,0.00001,0.09250,0.12500,0.00001,0.00001,0.12500,0.00001,0.12500,0.00001,0.00001,0.12500,0.09250,0.00001,0.00001,0.09250,0.00001,0.09250,0.00001,0.00001,0.09250,0.00003,0.09250,0.00001,0.00001,0.12500,0.00001,0.00001,0.09250,0.00001,0.00001,0.00001,0.12500,0.00001,0.00001,0.25000,0.00001,0.00001,0.12500,0.00001,0.00001,0.00001,0.09250,0.00001,0.00001,0.12500,0.00001,0.00001,0.09250,0.00000,0.09250,0.00001,0.00001,0.00001,0.12500,0.00001,0.00001,0.00001,0.09250,0.12500,0.00001,0.00001,0.00001,0.25000,0.00001,0.00001,0.00001,0.12500,0.09250,0.00001,0.00001,0.00001,0.12500,0.00001,0.00001,0.00001,0.09250,0.00000),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      states = c("normal", "iso_fat", "iso_mat", "het_fat", "het_mat"),
      symbols = c(
        "111", "112", "113", "121", "122", "123", "131", "132", "133",
        "211", "212", "213", "221", "222", "223", "231", "232", "233",
        "311", "312", "313", "321", "322", "323", "331", "332", "333",
        "000"
      )
    )
  )

  )
  
  



expected_def_blocks <- data.frame(
  seqnames = "6",
  start = 32489853,
  end = 32489876,
  group = "iso_mat",
  n_snps = 3,
  log_likelihood = 22.7,
  p_value = 1e-6,
  n_mendelian_error=2
)


file <- system.file(package = "UPDhmm", "extdata", "test.vcf.gz")
input <- VariantAnnotation::readVcf(file)
colnames(input) <- c("proband", "father", "mother")
test_that("Test if the general function works with a different hmm", {
  out <- calculateEvents(input,hmm = new_hmm)
  out$seqnames <- as.character(out$seqnames)
  out$log_likelihood <- floor(as.numeric(out$log_likelihood) * 10) / 10
  out$p_value <- floor(as.numeric(out$p_value) * 10^6) / 10^6
  out <- as.data.frame(out)
  expect_equal(expected_def_blocks, out)
  expect_s3_class(out, "data.frame")
})



  
  
  
  
  