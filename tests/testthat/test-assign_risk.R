ipssr <- function(karyotypes) {
  assign_risk(pk(karyotypes), scheme = "ipssr")
}

test_that("assign_risk: IPSS-R category table", {
  cases <- list(
    list("45,X,-Y", "Very Good"),
    list("46,XX,del(11)(q14q23)", "Very Good"),
    list("46,XX", "Good"),
    list("46,XX,del(5)(q13q33)", "Good"),
    list("46,XX,del(12)(p11p13)", "Good"),
    list("46,XX,del(20)(q11q13)", "Good"),
    list("47,XX,del(5)(q13q33),+8", "Good"),
    list("46,XX,del(7)(q22q36)", "Intermediate"),
    list("47,XX,+8", "Intermediate"),
    list("47,XX,+19", "Intermediate"),
    list("46,XX,i(17)(q10)", "Intermediate"),
    list("47,XX,+21", "Intermediate"),
    list("48,XX,+8,+21", "Intermediate"),
    list("45,XX,-7", "Poor"),
    list("46,XX,inv(3)(q21q26)", "Poor"),
    list("46,XX,t(3;3)(q21;q26)", "Poor"),
    list("46,XX,del(3)(q21q26)", "Poor"),
    list("46,XX,-7,+8", "Poor"),
    list("47,XX,del(7)(q22q36),+8", "Poor"),
    list("44,XX,-5,-7,del(5)(q13)", "Poor"),
    list("43,XX,-5,-7,-18,del(5)(q13)", "Very Poor")
  )
  for (tc in cases) {
    r <- ipssr(tc[[1]])
    expect_equal(as.character(r$ipssr_cyto_risk), tc[[2]], info = tc[[1]])
  }
})

test_that("assign_risk: score matches category", {
  r <- ipssr(c(
    "45,X,-Y",
    "46,XX",
    "47,XX,+8",
    "45,XX,-7",
    "43,XX,-5,-7,-18,del(5)(q13)"
  ))
  expect_equal(r$ipssr_cyto_score, c(0L, 1L, 2L, 3L, 4L))
})

test_that("assign_risk: risk column is an ordered factor", {
  r <- ipssr("46,XX")
  expect_s3_class(r$ipssr_cyto_risk, "ordered")
  expect_equal(
    levels(r$ipssr_cyto_risk),
    c("Very Good", "Good", "Intermediate", "Poor", "Very Poor")
  )
  expect_true(
    ipssr("45,XX,-7")$ipssr_cyto_risk > ipssr("46,XX")$ipssr_cyto_risk
  )
})

test_that("assign_risk: a double with both del(5q) and -7/del(7q) resolves to Poor", {
  expect_equal(
    as.character(ipssr("44,XX,del(5)(q13q33),-7")$ipssr_cyto_risk),
    "Poor"
  )
  expect_equal(
    as.character(ipssr("46,XX,del(5)(q13q33),del(7)(q22q36)")$ipssr_cyto_risk),
    "Poor"
  )
})

test_that("assign_risk: unscoreable and unparsed rows are NA", {
  r <- ipssr("45,X")
  expect_true(is.na(r$ipssr_cyto_risk))
  expect_true(is.na(r$ipssr_cyto_score))

  r2 <- suppressWarnings(ipssr(c("46,XX", "mos 47,XXY[10]")))
  expect_equal(as.character(r2$ipssr_cyto_risk[1]), "Good")
  expect_true(is.na(r2$ipssr_cyto_risk[2]))
})

test_that("assign_risk: appends to the parsed tibble without dropping columns", {
  parsed <- pk("46,XX,del(5)(q13q33)")
  r <- assign_risk(parsed)
  expect_true(all(names(parsed) %in% names(r)))
  expect_equal(nrow(r), nrow(parsed))
  expect_equal(
    attr(r, "karyoparser_version"),
    attr(parsed, "karyoparser_version")
  )
})

test_that("assign_risk: errors on non-data-frame input", {
  expect_error(assign_risk("46,XX"), "must be a data frame")
})

test_that("assign_risk: errors when required columns are absent", {
  narrowed <- pk("46,XX", columns = "general")
  expect_error(assign_risk(narrowed), "missing column\\(s\\) required")
  expect_error(assign_risk(narrowed), "distinct_aberrations")
})

test_that("assign_risk: rejects an unknown scheme", {
  expect_error(assign_risk(pk("46,XX"), scheme = "eln2022"))
})

test_that("assign_risk: empty input returns zero rows with both columns", {
  r <- assign_risk(pk(character(0)))
  expect_equal(nrow(r), 0L)
  expect_true(all(
    c("ipssr_cyto_risk", "ipssr_cyto_score") %in% names(r)
  ))
})
