test_that("example_karyotypes has expected structure", {
  expect_s3_class(example_karyotypes, "tbl_df")
  expect_named(example_karyotypes, c("sample_id", "karyotype"))
  expect_equal(nrow(example_karyotypes), 102L)
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
  expect_equal(nrow(result), 102L)
})

test_that("example_karyotypes covers every rule and general flag", {
  result <- suppressMessages(parse_karyo(
    example_karyotypes,
    karyotype_column = "karyotype",
    id_column = "sample_id",
    on_issues = "preprocess",
    verbose = FALSE
  ))
  # Keyed on the catalog, not on myeloid_rules alone, so a newly added rule or
  # general flag is guarded without remembering to widen this test.
  catalog <- karyoparser:::.column_catalog(myeloid_rules)
  flags <- catalog$column[catalog$class %in% c("rule", "general")]
  uncovered <- flags[vapply(
    flags,
    \(f) sum(result[[f]], na.rm = TRUE) == 0L,
    logical(1)
  )]
  expect_equal(uncovered, character(0))
})
