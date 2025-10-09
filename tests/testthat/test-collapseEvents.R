# Input test dataframe
test_df <- data.frame(
    ID = c("S1", "S1", "S1", "S2", "S2"),
    seqnames = c("1", "1", "1", "2", "2"),
    start = c(100, 150, 300, 500, 550),
    end = c(120, 180, 320, 520, 580),
    group = c("iso_mat", "iso_mat", "het_pat", "iso_mat", "iso_mat"),
    n_mendelian_error = c(5, 10, 2, 50, 30),
    stringsAsFactors = FALSE
  )
  
  # Expected result after collapsing
expected_result <- data.frame(
    ID = c("S1", "S1", "S2"),
    seqnames = c("1", "1", "2"),
    group = c("het_pat","iso_mat", "iso_mat"),
    n_events = c(1, 2, 2),
    total_mendelian_error = c(2, 15, 80),
    total_size = c(20, 50, 50),
    collapsed_events = c( "1:300-320","1:100-120,1:150-180", "2:500-520,2:550-580"),
    min_start = c(300,100, 500),
    max_end = c(320, 180, 580),
    stringsAsFactors = FALSE
  )
  



test_that("Test if calculation collapseEvents works correctly", {
  # Run function
  out <- collapseEvents(test_df)
  # Test equality
  expect_equal(out, expected_result)
})



