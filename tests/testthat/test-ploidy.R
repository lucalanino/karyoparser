test_that("ploidy_from_count covers all thresholds", {
  expect_equal(ploidy_from_count(25), "near_haploid")
  expect_equal(ploidy_from_count(23), "near_haploid")
  expect_equal(ploidy_from_count(29), "near_haploid")
  expect_equal(ploidy_from_count(30), "low_hypodiploid")
  expect_equal(ploidy_from_count(33), "low_hypodiploid")
  expect_equal(ploidy_from_count(35), "other")
  expect_equal(ploidy_from_count(40), "high_hypodiploid")
  expect_equal(ploidy_from_count(45), "high_hypodiploid")
  expect_equal(ploidy_from_count(46), "diploid")
  expect_equal(ploidy_from_count(48), "other")
  expect_equal(ploidy_from_count(51), "hyperdiploid")
  expect_equal(ploidy_from_count(92), "hyperdiploid")
  expect_equal(ploidy_from_count(NA), "unknown")
})

test_that("ploidy_category for composite karyotype picks most abnormal", {
  res <- ploidy_category("46,XX[10]/42,XX,-5,-7,-8,-9[10]")
  expect_equal(res$ploidy, "high_hypodiploid")
  expect_equal(res$chromosome_count, 42L)
})

test_that("ploidy_category filters clones with <5 metaphases", {
  res <- ploidy_category("25,X[3]/46,XX[10]")
  expect_equal(res$ploidy, "diploid")
  expect_equal(res$chromosome_count, 46L)
})

test_that("ploidy_category detects mixed ploidy", {
  res <- ploidy_category("45,XX,-7[10]/51,XX,+8,+11,+13,+19,+21[8]")
  expect_true(res$mixed)
})

test_that("ploidy_category no brackets uses all clones", {
  res <- ploidy_category("46,XX/45,XY,-7")
  expect_equal(res$ploidy, "high_hypodiploid")
  expect_equal(res$chromosome_count, 45L)
})

test_that("ploidy_category: range metaphase bracket [N~M] uses max of range for eligibility", {
  res <- ploidy_category("25,X[3]/47,XX,+8[3~5]")
  expect_equal(res$ploidy, "other")
  expect_equal(res$chromosome_count, 47L)
})

test_that("ploidy_from_count: values below near_haploid floor return 'other'", {
  expect_equal(ploidy_from_count(22), "other")
})

test_that("ploidy_from_count: exact gap boundaries return 'other'", {
  expect_equal(ploidy_from_count(34), "other")
  expect_equal(ploidy_from_count(39), "other")
  expect_equal(ploidy_from_count(47), "other")
  expect_equal(ploidy_from_count(50), "other")
})

test_that("ploidy_category: all clones with <5 metaphases returns 'unknown'", {
  res <- ploidy_category("46,XX[2]/45,XY,-7[3]")
  expect_equal(res$ploidy, "unknown")
  expect_true(is.na(res$chromosome_count))
  expect_false(res$mixed)
})

test_that("ploidy_category: range chromosome count uses round(mean())", {
  res <- ploidy_category("45~46,XX")
  expect_equal(res$ploidy, "diploid")
  expect_equal(res$chromosome_count, 46L)
})
