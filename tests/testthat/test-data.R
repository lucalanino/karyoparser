test_that("example_karyotypes has expected structure", {
  expect_s3_class(example_karyotypes, "tbl_df")
  expect_named(example_karyotypes, c("sample_id", "karyotype"))
  expect_equal(nrow(example_karyotypes), 100L)
  expect_false(anyDuplicated(example_karyotypes$sample_id) > 0L)
})

test_that("example_karyotypes parses end-to-end without error", {
  result <- parse_karyo(
    example_karyotypes,
    karyotype_column = "karyotype",
    id_column = "sample_id",
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_equal(nrow(result), 100L)
})
