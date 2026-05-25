test_that("pk() is available and returns a tibble", {
  r <- pk("46,XX")
  expect_s3_class(r, "tbl_df")
})

test_that("pk() passes on_issues = 'warn': warns on issues, does not error", {
  expect_warning(
    r <- pk(c("46,XX", "not_a_karyotype")),
    regexp = "had issues"
  )
  expect_equal(nrow(r), 2L)
})

test_that("pk() passes verbose = FALSE: no messages emitted", {
  expect_no_message(pk("46,XX"))
})

test_that("pk() forwards extra args to parse_karyo()", {
  r <- pk("46,XX", rules = myeloid_rules)
  expect_s3_class(r, "tbl_df")
})
