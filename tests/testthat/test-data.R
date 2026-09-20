test_that("example_karyotypes has expected structure", {
  expect_s3_class(example_karyotypes, "tbl_df")
  expect_named(example_karyotypes, c("sample_id", "karyotype"))
  expect_equal(nrow(example_karyotypes), 101L)
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
  expect_equal(nrow(result), 101L)
})

test_that("example_karyotypes covers every myeloid_rules flag", {
  result <- suppressMessages(parse_karyo(
    example_karyotypes,
    karyotype_column = "karyotype",
    id_column = "sample_id",
    on_issues = "preprocess",
    verbose = FALSE
  ))
  flags <- karyoparser:::.sanitize_flag_name(myeloid_rules$flag_name)
  uncovered <- flags[vapply(
    flags,
    \(f) sum(result[[f]], na.rm = TRUE) == 0L,
    logical(1)
  )]
  expect_equal(uncovered, character(0))
})
