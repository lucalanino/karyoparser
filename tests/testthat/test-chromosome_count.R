test_that("chromosome_count_from_karyotype for composite picks most abnormal", {
  expect_equal(
    chromosome_count_from_karyotype("46,XX[10]/42,XX,-5,-7,-8,-9[10]"),
    42L
  )
})

test_that("chromosome_count_from_karyotype default min_metaphases is 2", {
  expect_equal(chromosome_count_from_karyotype("25,X[3]/46,XX[10]"), 25L)
  expect_equal(chromosome_count_from_karyotype("25,X[1]/46,XX[10]"), 46L)
})

test_that("chromosome_count_from_karyotype filters clones below min_metaphases", {
  expect_equal(
    chromosome_count_from_karyotype("25,X[3]/46,XX[10]", min_metaphases = 5),
    46L
  )
})

test_that("chromosome_count_from_karyotype min_metaphases = 0 keeps every counted clone", {
  expect_equal(
    chromosome_count_from_karyotype("25,X[1]/46,XX[10]", min_metaphases = 0),
    25L
  )
})

test_that("chromosome_count_from_karyotype no brackets uses all clones", {
  expect_equal(chromosome_count_from_karyotype("46,XX/45,XY,-7"), 45L)
})

test_that("chromosome_count_from_karyotype no brackets ignores min_metaphases", {
  expect_equal(
    chromosome_count_from_karyotype("46,XX/45,XY,-7", min_metaphases = 50),
    45L
  )
})

test_that("chromosome_count_from_karyotype: [N~M] uses max of range for eligibility", {
  expect_equal(
    chromosome_count_from_karyotype(
      "25,X[3]/47,XX,+8[3~5]",
      min_metaphases = 5
    ),
    47L
  )
})

test_that("chromosome_count_from_karyotype: all clones below threshold returns NA", {
  expect_true(is.na(chromosome_count_from_karyotype("46,XX[1]/45,XY,-7[1]")))
  expect_true(is.na(
    chromosome_count_from_karyotype("46,XX[2]/45,XY,-7[3]", min_metaphases = 5)
  ))
})

test_that("chromosome_count_from_karyotype: range count uses round(mean())", {
  expect_equal(chromosome_count_from_karyotype("45~46,XX"), 46L)
})
