test_that("unicode normalization: NBSPs detected and fixed by preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo("46, XX"))
  expect_equal(result$preprocessed, "46,XX")
  expect_equal(result$status, "fixed")
  r_warn <- suppressWarnings(pk("46, XX"))
  expect_true(is.na(r_warn$normal_karyotype))
  r_fix <- parse_karyo("46, XX", on_issues = "preprocess", verbose = FALSE)
  expect_equal(r_fix$normal_karyotype, 1L)
})

test_that("whitespace and delimiter cleanup: normalize_iscn runs inside preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo("46 , XX , +8"))
  expect_equal(result$preprocessed, "46,XX,+8")
  r <- parse_karyo("46 , XX , +8", on_issues = "preprocess", verbose = FALSE)
  expect_equal(r$tris8, 1L)
})

test_that("psu dic normalization: preprocess_karyo normalizes before parse", {
  result <- suppressMessages(preprocess_karyo("46,XX,PSU DIC(15;22)"))
  expect_equal(result$preprocessed, "46,XX,psu dic(15;22)")
  r <- parse_karyo(
    "46,XX,PSU DIC(15;22)",
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_equal(r$pseudodicentric, 1L)
})

test_that("idem and sl case normalization via preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo(
    "46,XX,del(5)(q13)[10]/46,IDEM,+8[5]"
  ))
  expect_equal(result$preprocessed, "46,XX,del(5)(q13)[10]/46,idem,+8[5]")
  r <- parse_karyo(
    "46,XX,del(5)(q13)[10]/46,IDEM,+8[5]",
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_equal(r$comma_count_aberrations, 2L)
})

test_that("duplicate commas and leading/trailing delimiters: normalize_iscn in preprocess", {
  result <- suppressMessages(preprocess_karyo(",46,,XX,"))
  expect_equal(result$preprocessed, "46,XX")
  r <- parse_karyo(",46,,XX,", on_issues = "preprocess", verbose = FALSE)
  expect_equal(r$normal_karyotype, 1L)
})

test_that("fullwidth plus: detected as unicode_notation, fixed by preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo("47,XY,＋8"))
  expect_equal(result$preprocessed, "47,XY,+8")
  r <- parse_karyo("47,XY,＋8", on_issues = "preprocess", verbose = FALSE)
  expect_equal(r$tris8, 1L)
})

test_that("cp bracket spacing normalization via parse_karyo", {
  r <- pk("46,XX[cp 20]")
  expect_true(!is.na(r$total_metaphases))
})

test_that("range count: dash normalized to tilde via parse_karyo", {
  r <- pk("46-47,XY,+8")
  expect_equal(r$normal_karyotype, 0L)
})

test_that("karyotype starting with 45 correctly assigns monosomy", {
  r <- pk("45,XY,-7")
  expect_equal(r$mono7, 1L)
})

test_that("preprocess_karyo: returns tibble with original/preprocessed/status columns", {
  result <- suppressMessages(preprocess_karyo(
    "46,XX,t(9;22)(q34;q11)[20] .Clinical note here"
  ))
  expect_true(tibble::is_tibble(result))
  expect_named(result, c("original", "preprocessed", "status"))
})

test_that("preprocess_karyo: strips trailing narrative after bracket", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,t(9;22)(q34;q11)[20] .Clinical note here"
    ))$preprocessed,
    "46,XX,t(9;22)(q34;q11)[20]"
  )
})

test_that("preprocess_karyo: strips trailing narrative after bracket with no space", {
  expect_equal(
    suppressMessages(preprocess_karyo("47,XY,+21[10].Some text"))$preprocessed,
    "47,XY,+21[10]"
  )
})

test_that("preprocess_karyo: decodes HTML lt/gt entities", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,t(9;22)(q34;q11) &lt;AML&gt;"
    ))$preprocessed,
    "46,XX,t(9;22)(q34;q11)<AML>"
  )
})

test_that("preprocess_karyo: &amp; entity decoded but decoded form is unfixable (no_sex_complement remains)", {
  result <- suppressMessages(preprocess_karyo("46,XX &amp; notes"))
  expect_equal(result$status, "unfixable")
  expect_equal(result$preprocessed, NA_character_)
})

test_that("preprocess_karyo: .// and // prefix are unfixable; preprocessed is NA", {
  result <- suppressMessages(preprocess_karyo(c(".//46,XX", "..//46,XX")))
  expect_equal(result$preprocessed, c(NA_character_, NA_character_))
  expect_equal(result$status, c("unfixable", "unfixable"))
})

test_that("preprocess_karyo: strips leading dot before digit", {
  result <- suppressMessages(preprocess_karyo(c(".46,XX", "..47,XY,+21")))
  expect_equal(result$preprocessed, c("46,XX", "47,XY,+21"))
  expect_equal(result$status, c("fixed", "fixed"))
})

test_that("preprocess_karyo: leaves clean strings unchanged", {
  clean <- c("46,XX", "47,XY,+21[10]", "46,XX,t(9;22)(q34;q11)[20]")
  result <- suppressMessages(preprocess_karyo(clean))
  expect_equal(result$preprocessed, clean)
  expect_equal(result$status, rep("clean", 3L))
})

test_that("preprocess_karyo: idempotent on preprocessed column", {
  dirty <- c(".46,XX", ".//47,XY,+21", "46,XX[10] .Note")
  once <- suppressMessages(preprocess_karyo(dirty))
  twice <- suppressMessages(preprocess_karyo(once$preprocessed))
  expect_equal(once$preprocessed, twice$preprocessed)
})

test_that("preprocess_karyo: does not corrupt mid-string .add() tokens", {
  x <- "46,XX,.add(1)(q21)"
  expect_equal(suppressMessages(preprocess_karyo(x))$preprocessed, x)
})

test_that("preprocess_karyo: handles empty vector", {
  result <- preprocess_karyo(character(0))
  expect_true(tibble::is_tibble(result))
  expect_equal(nrow(result), 0L)
  expect_named(result, c("original", "preprocessed", "status"))
})

test_that("fish_notation: detected when nuc ish suffix present", {
  result <- suppressMessages(check_karyo(
    "[12]/46,XX[3] .nuc ish(PDGFRA x3)[20/200]"
  ))
  expect_equal(result$fish_notation, 1L)
})

test_that("fish_notation: detected when .ish suffix present", {
  result <- suppressMessages(check_karyo("46,XX,t(9;22)[15] .ish(BCR-ABL)"))
  expect_equal(result$fish_notation, 1L)
})

test_that("fish_notation: not detected for clean karyotype", {
  result <- suppressMessages(check_karyo("46,XX,t(9;22)(q34;q11)[15]/46,XX[5]"))
  expect_equal(result$fish_notation, 0L)
})

test_that("preprocess_karyo: strips nuc ish suffix after bracket", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,t(9;22)[15]/46,XX[3] .nuc ish(PDGFRA x3)[20/200]"
    ))$preprocessed,
    "46,XX,t(9;22)[15]/46,XX[3]"
  )
})

test_that("preprocess_karyo: strips .ish suffix after bracket", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,t(9;22)(q34;q11)[15] .ish(BCR-ABL)"
    ))$preprocessed,
    "46,XX,t(9;22)(q34;q11)[15]"
  )
})

test_that("fish_notation: detected for ).ish attachment", {
  result <- suppressMessages(check_karyo(
    "46,XX,der(17)t(11;17)(q14;q11).ish der(15)t(15;17)(RARA+,PML+)"
  ))
  expect_equal(result$fish_notation, 1L)
})

test_that("fish_notation: detected for dmin.ish attachment", {
  result <- suppressMessages(check_karyo(
    "46,XY,10~>50dmin.ish del(8)(q24q24)(MYC-),dmin(MYC+)"
  ))
  expect_equal(result$fish_notation, 1L)
})

test_that("preprocess_karyo: strips .ish suffix after closing paren", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,der(17)t(11;17)(q14;q11).ish der(15)t(15;17)(RARA+,PML+),der(17)t(11;17)(RARA-,PML-)"
    ))$preprocessed,
    "46,XX,der(17)t(11;17)(q14;q11)"
  )
})

test_that("preprocess_karyo: strips dmin.ish suffix, preserves dmin count", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XY,10~>50dmin.ish del(8)(q24q24)(MYC-),dmin(MYC+)"
    ))$preprocessed,
    "46,XY,10~>50dmin"
  )
})

test_that("preprocess_karyo: strips cp.ish suffix, preserves cp notation", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "40~43,XX,der(21)t(19;21)(p12;q22)x1~3 cp.ish der(21)t(19;21)(RUNX1+)"
    ))$preprocessed,
    "40~43,XX,der(21)t(19;21)(p12;q22)x1~3 cp"
  )
})

test_that("preprocess_karyo: strips .ish after aberration, preserves aberration", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,+22.ish der(3;21)(5'MECOM+),der(18)(5'MECOM+)"
    ))$preprocessed,
    "46,XX,+22"
  )
})

test_that("mar_space: detected when space between count and mar", {
  result <- suppressMessages(check_karyo("47,XY,+1~4 mar[cp15]"))
  expect_equal(result$mar_space, 1L)
})

test_that("mar_space: not detected for well-formed mar token", {
  result <- suppressMessages(check_karyo("47,XY,+mar[5]"))
  expect_equal(result$mar_space, 0L)
})

test_that("preprocess_karyo: removes space before mar token", {
  expect_equal(
    suppressMessages(preprocess_karyo("47,XY,+1~4 mar[cp15]"))$preprocessed,
    "47,XY,+1~4mar[cp15]"
  )
})

test_that("preprocess_karyo: mar_space fix handles minus count", {
  expect_equal(
    suppressMessages(preprocess_karyo("46,XX,-1 mar[10]"))$preprocessed,
    "46,XX,-1mar[10]"
  )
})

test_that("midstring_linewrap: detected for ', .+N' artifact", {
  result <- suppressMessages(check_karyo("46,XY,del(5)(q13), .+8[10]/46,XY[5]"))
  expect_equal(result$midstring_linewrap, 1L)
})

test_that("preprocess_karyo: collapses ', .+8' mid-string artifact", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XY,del(5)(q13), .+8[10]/46,XY[5]"
    ))$preprocessed,
    "46,XY,del(5)(q13),+8[10]/46,XY[5]"
  )
})

test_that("midstring_linewrap: not triggered by trailing narrative with uppercase", {
  expect_equal(
    suppressMessages(preprocess_karyo("46,XX[20] .Abnormal note"))$preprocessed,
    "46,XX[20]"
  )
})
