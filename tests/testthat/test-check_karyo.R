test_that("check_karyo: clean input returns one row per string, all zeros", {
  result <- suppressMessages(check_karyo(c("46,XX", "47,XY,+21[10]")))
  expect_equal(nrow(result), 2L)
  expect_equal(result$fixable, c(0L, 0L))
  expect_equal(result$unfixable, c(0L, 0L))
})

test_that("check_karyo: column schema -- fixable/unfixable first, then sorted unfixable issue cols", {
  result <- suppressMessages(check_karyo("46,XX"))
  unfixable_cols <- sort(setdiff(
    karyoparser:::.all_issue_types,
    karyoparser:::.fixable_issue_types
  ))
  expected_cols <- c(
    "karyotype",
    "fixable",
    "unfixable",
    unfixable_cols
  )
  expect_named(result, expected_cols)
})

test_that("check_karyo: empty/NA detected as unfixable", {
  result <- suppressMessages(check_karyo(c(NA, "")))
  expect_equal(result$empty, c(1L, 1L))
  expect_equal(result$unfixable, c(1L, 1L))
})

test_that("check_karyo: no chromosome count detected", {
  result <- suppressMessages(check_karyo("XX,+8"))
  expect_equal(result$no_chromosome_count, 1L)
  expect_equal(result$unfixable, 1L)
})

test_that("check_karyo: unbalanced parentheses detected", {
  result <- suppressMessages(check_karyo("46,XX,del(5)(q13"))
  expect_equal(result$unbalanced_parentheses, 1L)
  expect_equal(result$unfixable, 1L)
})

test_that("check_karyo: unbalanced brackets detected", {
  result <- suppressMessages(check_karyo("46,XX[10"))
  expect_equal(result$unbalanced_brackets, 1L)
  expect_equal(result$unfixable, 1L)
})

test_that("check_karyo: unbalanced brackets in a later clone survive trailing_narrative's fix", {
  result <- suppressMessages(check_karyo("46,XX[10]/46,XX,del(5)(q13)[5"))
  expect_equal(result$unbalanced_brackets, 1L)
  expect_equal(result$unfixable, 1L)
})

test_that("check_karyo: dirty markers detected as fixable", {
  result <- suppressMessages(check_karyo(".47,XY,+21"))
  expect_equal(has_issue(".47,XY,+21", "leading_dot"), 1L)
  expect_equal(result$fixable, 1L)
  expect_equal(result$unfixable, 0L)
})

test_that("check_karyo: trailing narrative detected as fixable", {
  result <- suppressMessages(check_karyo("46,XX[20] .some text"))
  expect_equal(has_issue("46,XX[20] .some text", "trailing_narrative"), 1L)
  expect_equal(result$fixable, 1L)
})

test_that("check_karyo: html entities detected as fixable", {
  k <- "46,XX,t(9;22)(q34;q11) &lt;AML&gt;"
  result <- suppressMessages(check_karyo(k))
  expect_equal(has_issue(k, "html_entities"), 1L)
  expect_equal(result$fixable, 1L)
})

test_that("check_karyo: Greek homoglyph sex complement detected as fixable, not no_sex_complement", {
  k <- "45,\u03A7\u03A7,-7[22]/47,\u03A7\u03A7,+8[3]"
  result <- suppressMessages(check_karyo(k))
  expect_equal(has_issue(k, "non_ascii_homoglyph"), 1L)
  expect_equal(result$no_sex_complement, 0L)
  expect_equal(result$fixable, 1L)
  expect_equal(result$unfixable, 0L)
})

test_that("check_karyo: stray non-ASCII character with no known fix is unfixable", {
  result <- suppressMessages(check_karyo(
    "46,XX,der(1)t(1;7)(p11;p11)\u00B5"
  ))
  expect_equal(result$stray_non_ascii, 1L)
  expect_equal(result$unfixable, 1L)
})

test_that("check_karyo: one row per input; clean/dirty rows correctly flagged", {
  x <- c("46,XX", ".47,XY,+21", "46,XX[5] .note", "46,XY")
  result <- suppressMessages(check_karyo(x))
  expect_equal(nrow(result), 4L)
  expect_equal(result$fixable, c(0L, 1L, 1L, 0L))
})

test_that("check_karyo: multiple issues on same row all flagged", {
  k <- ".46,XX[10] .note"
  result <- suppressMessages(check_karyo(k))
  expect_equal(has_issue(k, "leading_dot"), 1L)
  expect_equal(has_issue(k, "trailing_narrative"), 1L)
  expect_equal(result$fixable, 1L)
})

test_that("check_karyo: NA input detected as unfixable empty", {
  result <- suppressMessages(check_karyo(c("46,XX", NA, "47,XY,+21")))
  expect_equal(nrow(result), 3L)
  expect_equal(result$empty[2], 1L)
  expect_equal(result$unfixable[2], 1L)
  expect_equal(result$fixable[1], 0L)
  expect_equal(result$fixable[3], 0L)
})

test_that("check_karyo: empty vector returns zero-row tibble with full schema", {
  result <- check_karyo(character(0))
  expect_equal(nrow(result), 0L)
  expect_true("fixable" %in% names(result))
  expect_true("unfixable" %in% names(result))
})

test_that("check_karyo detects invalid_idem when idem appears in clone 1", {
  result <- suppressMessages(check_karyo("46,XX,idem[10]"))
  expect_equal(result$invalid_idem, 1L)
  expect_equal(result$unfixable, 1L)
})

test_that("check_karyo detects unparseable_bracket for non-numeric bracket content", {
  result <- suppressMessages(check_karyo("46,XX[abc]"))
  expect_equal(result$unparseable_bracket, 1L)
  expect_equal(result$unfixable, 1L)
})

test_that("XXYY and XXXY are valid sex complements (no no_sex_complement fired)", {
  expect_equal(
    suppressMessages(check_karyo("48,XXYY,+1[10]"))$no_sex_complement,
    0L
  )
  expect_equal(
    suppressMessages(check_karyo("48,XXXY,+1[10]"))$no_sex_complement,
    0L
  )
})

test_that("parse_karyo handles 48,XXYY karyotype", {
  r <- pk("48,XXYY,+1[10]")
  expect_false(is.na(r$chromosome_count))
  expect_equal(r$tris1, 1L)
})

test_that("check_karyo: no_sex_complement detected when sex chromosome token absent", {
  result <- suppressMessages(check_karyo("46,+8"))
  expect_equal(result$no_sex_complement, 1L)
  expect_equal(result$unfixable, 1L)
})

test_that("no_sex_complement: clone with only sex-chromosome aberrations is accepted", {
  for (k in c(
    "51,add(X)(q26),-Y,+5,+7[10]",
    "44,-X,t(X;14)(q28;q11.2),-2,+mar[10]",
    "45,-Xx2,t(7;15)(p13;q11.2),-9[10]",
    "41~43,del(X)(q24),der(X)t(X;11)(q22;q13),-12[10]"
  )) {
    chk <- suppressMessages(check_karyo(k))
    expect_equal(chk$no_sex_complement, 0L, info = k)
    expect_equal(chk$unfixable, 0L, info = k)
  }
})

test_that("mosaic_karyotype flagged for 'mos' prefix, not no_chromosome_count", {
  result <- suppressMessages(check_karyo("mos 47,XXY[10]/46,XY[5]"))
  expect_equal(result$mosaic_karyotype, 1L)
  expect_equal(result$no_chromosome_count, 0L)
  expect_equal(result$unfixable, 1L)
})

test_that("non_clonal_sca flagged for ncSCA token (leading or mid-clone), fixable", {
  lead_k <- "ncSCA[4]/46,XY[11]"
  lead <- suppressMessages(check_karyo(lead_k))
  expect_equal(has_issue(lead_k, "non_clonal_sca"), 1L)
  expect_equal(lead$no_chromosome_count, 0L)
  expect_equal(lead$fixable, 1L)
  expect_equal(lead$unfixable, 0L)

  mid_k <- "46,XX(ncSCA)[1]//46,XY[19]"
  mid <- suppressMessages(check_karyo(mid_k))
  expect_equal(has_issue(mid_k, "non_clonal_sca"), 1L)
  expect_equal(mid$no_sex_complement, 0L)
  expect_equal(mid$fixable, 1L)
  expect_equal(mid$unfixable, 0L)
})

test_that("no_chromosome_count still fires for genuine non-count strings", {
  for (k in c("XX,+8", "Abnormal clones detected", "FALSE")) {
    r <- suppressMessages(check_karyo(k))
    expect_equal(r$no_chromosome_count, 1L, info = k)
    expect_equal(r$mosaic_karyotype, 0L, info = k)
    expect_equal(has_issue(k, "non_clonal_sca"), 0L, info = k)
  }
})

test_that("no_sex_complement: bare count+aberration with no sex reference still fires", {
  # '+8' carries no sex-chromosome operator, so the gate must hold.
  expect_equal(suppressMessages(check_karyo("46,+8"))$no_sex_complement, 1L)
  expect_equal(suppressMessages(check_karyo("46,+8"))$unfixable, 1L)
})

test_that("single_token fires for a bare comma-less string", {
  # E.g. Excel silently drops the '+' from a numeric-looking '+8' cell,
  # leaving '8' indistinguishable from a bare chromosome count.
  for (k in c("8", "46", "22")) {
    result <- suppressMessages(check_karyo(k))
    expect_equal(result$single_token, 1L, info = k)
    expect_equal(result$unfixable, 1L, info = k)
    expect_equal(result$no_chromosome_count, 0L, info = k)
    expect_equal(result$no_sex_complement, 0L, info = k)
  }
})

test_that("single_token does not fire once a comma is present", {
  for (k in c("46,XX", "46,+8", "8,del(5)(q13q33)")) {
    result <- suppressMessages(check_karyo(k))
    expect_equal(result$single_token, 0L, info = k)
  }
})

test_that("single_token yields to no_chromosome_count for signed bare tokens", {
  # '-8'/'+8' don't start with a digit, so they're still caught earlier.
  for (k in c("-8", "+8")) {
    result <- suppressMessages(check_karyo(k))
    expect_equal(result$no_chromosome_count, 1L, info = k)
    expect_equal(result$single_token, 0L, info = k)
    expect_equal(result$unfixable, 1L, info = k)
  }
})

test_that("constitutional sex complement is flagged unfixable, not no_sex_complement", {
  for (k in c(
    "47,XXYc[20]",
    "47,XXXc[20]",
    "47,XXYc?[20]",
    "47,XXYc[14]/46,XX,-Y[6]"
  )) {
    result <- suppressMessages(check_karyo(k))
    expect_equal(result$constitutional_sex_complement, 1L, info = k)
    expect_equal(result$no_sex_complement, 0L, info = k)
    expect_equal(result$unfixable, 1L, info = k)
  }
  # A plain (non-constitutional) complement is unaffected.
  ok <- suppressMessages(check_karyo("47,XXY[20]"))
  expect_equal(ok$constitutional_sex_complement, 0L)
  expect_equal(ok$unfixable, 0L)
})
