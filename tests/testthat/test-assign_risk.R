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
  expect_error(assign_risk(pk("46,XX"), scheme = "ipssm"))
})

test_that("assign_risk: empty input returns zero rows with both columns", {
  r <- assign_risk(pk(character(0)))
  expect_equal(nrow(r), 0L)
  expect_true(all(
    c("ipssr_cyto_risk", "ipssr_cyto_score") %in% names(r)
  ))
})

eln <- function(karyotypes) {
  assign_risk(pk(karyotypes), scheme = "eln2022")
}

test_that("assign_risk: ELN 2022 category table", {
  cases <- list(
    list("46,XX,t(8;21)(q22;q22)", "Favorable"),
    list("46,XX,inv(16)(p13q22)", "Favorable"),
    list("46,XX,t(16;16)(p13;q22)", "Favorable"),
    list("46,XX", "Intermediate"),
    list("47,XX,+8", "Intermediate"),
    list("46,XX,del(20)(q11q13)", "Intermediate"),
    list("46,XX,t(6;9)(p22;q34)", "Adverse"),
    list("46,XX,t(9;22)(q34;q11)", "Adverse"),
    list("46,XX,t(8;16)(p11;p13)", "Adverse"),
    list("46,XX,inv(3)(q21q26)", "Adverse"),
    list("46,XX,t(3;3)(q21;q26)", "Adverse"),
    list("46,XX,t(3;21)(q26;q22)", "Adverse"),
    list("45,XX,-5", "Adverse"),
    list("46,XX,del(5)(q13q33)", "Adverse"),
    list("45,XX,-7", "Adverse"),
    list("45,XX,-17", "Adverse"),
    list("46,XX,del(17)(p11p13)", "Adverse"),
    list("46,XX,add(17)(p11)", "Adverse"),
    list("46,XX,i(17)(q10)", "Adverse"),
    list("46,XX,del(20)(q11),del(13)(q14),+8", "Adverse")
  )
  for (tc in cases) {
    r <- eln(tc[[1]])
    expect_equal(as.character(r$eln2022_cyto_risk), tc[[2]], info = tc[[1]])
  }
})

test_that("assign_risk: ELN t(6;9) is adverse at the ISCN p23 spelling too", {
  expect_equal(
    as.character(eln("46,XX,t(6;9)(p23;q34)")$eln2022_cyto_risk),
    "Adverse"
  )
})

test_that("assign_risk: ELN t(9;11) carve-out beats the adverse KMT2A row", {
  expect_equal(pk("46,XX,t(9;11)(p21;q23)")$t_v_11q23, 1L)
  expect_equal(
    as.character(eln("46,XX,t(9;11)(p21;q23)")$eln2022_cyto_risk),
    "Intermediate"
  )
  expect_equal(
    as.character(eln("45,XX,-7,t(9;11)(p21;q23)")$eln2022_cyto_risk),
    "Intermediate"
  )
  expect_equal(
    as.character(eln("46,XX,t(11;19)(q23;p13)")$eln2022_cyto_risk),
    "Adverse"
  )
})

test_that("assign_risk: ELN favorable lesions outrank adverse criteria", {
  expect_equal(
    as.character(eln("43,XX,-5,-7,-18,t(8;21)(q22;q22)")$eln2022_cyto_risk),
    "Favorable"
  )
})

test_that("assign_risk: ELN hyperdiploid carve-out excludes pure-gain complex", {
  expect_equal(
    as.character(eln("49,XX,+8,+13,+21")$eln2022_cyto_risk),
    "Intermediate"
  )
  expect_equal(
    as.character(eln("52,XX,+1,+6,+8,+9,+14,+21")$eln2022_cyto_risk),
    "Intermediate"
  )
  expect_equal(
    as.character(eln("49,XX,+8,+13,+21,del(20)(q11)")$eln2022_cyto_risk),
    "Adverse"
  )
  expect_equal(
    as.character(eln("48,XX,+8,+13,-18")$eln2022_cyto_risk),
    "Adverse"
  )
  expect_equal(
    as.character(eln("49,XX,+8,+13,+21,+mar")$eln2022_cyto_risk),
    "Adverse"
  )
})

test_that("assign_risk: the carve-out does not leak into complex_karyotype", {
  r <- eln("49,XX,+8,+13,+21")
  expect_equal(r$complex_karyotype, 1L)
  expect_equal(as.character(r$eln2022_cyto_risk), "Intermediate")
})

test_that("assign_risk: eln2022 risk is an ordered factor and adds no score", {
  r <- eln("46,XX")
  expect_s3_class(r$eln2022_cyto_risk, "ordered")
  expect_equal(
    levels(r$eln2022_cyto_risk),
    c("Favorable", "Intermediate", "Adverse")
  )
  expect_false("eln2022_cyto_score" %in% names(r))
  expect_true(
    eln("45,XX,-7")$eln2022_cyto_risk > eln("46,XX")$eln2022_cyto_risk
  )
})

test_that("assign_risk: eln2022 returns NA for unparsed rows", {
  r <- suppressWarnings(eln(c("46,XX", "mos 47,XXY[10]")))
  expect_equal(as.character(r$eln2022_cyto_risk[1]), "Intermediate")
  expect_true(is.na(r$eln2022_cyto_risk[2]))
})

test_that("assign_risk: schemes chain and coexist", {
  r <- assign_risk(pk("45,XX,-7"), "ipssr") |> assign_risk("eln2022")
  expect_equal(as.character(r$ipssr_cyto_risk), "Poor")
  expect_equal(as.character(r$eln2022_cyto_risk), "Adverse")
})

test_that("assign_risk: eln2022 errors when required columns are absent", {
  narrowed <- pk("46,XX", columns = "classification")
  expect_error(
    assign_risk(narrowed, "eln2022"),
    "missing column\\(s\\) required"
  )
})
