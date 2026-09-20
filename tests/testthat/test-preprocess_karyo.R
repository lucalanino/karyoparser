test_that("unicode normalization: NBSPs detected and fixed by preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo("46,\u00A0XX"))
  expect_equal(result$preprocessed, "46,XX")
  expect_equal(result$status, "fixed")
  r_warn <- suppressWarnings(pk("46,\u00A0XX"))
  expect_true(is.na(r_warn$normal_karyotype))
  r_fix <- parse_karyo("46,\u00A0XX", on_issues = "preprocess", verbose = FALSE)
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
  expect_equal(r$general_pseudodicentric, 1L)
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
  result <- suppressMessages(preprocess_karyo("47,XY,\uFF0B8"))
  expect_equal(result$preprocessed, "47,XY,+8")
  r <- parse_karyo("47,XY,\uFF0B8", on_issues = "preprocess", verbose = FALSE)
  expect_equal(r$tris8, 1L)
})

test_that("zero-width spaces: stripped by preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo(
    "46,XX,t(3;21)(q26;q22)\u200B"
  ))
  expect_equal(result$preprocessed, "46,XX,t(3;21)(q26;q22)")
  expect_equal(result$status, "fixed")
})

test_that("Greek homoglyph sex complement: normalized to Latin X/Y by preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo(
    "45,\u03A7\u03A7,-7[22]/47,\u03A7\u03A7,+8[3]"
  ))
  expect_equal(result$preprocessed, "45,XX,-7[22]/47,XX,+8[3]")
  expect_equal(result$status, "fixed")
})

test_that("fullwidth ISCN punctuation: brackets/parens/comma normalized by preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo(
    "46\uFF0CXX\uFF0Cdel\uFF085\uFF09\uFF08q13\uFF09\uFF3B10\uFF3D"
  ))
  expect_equal(result$preprocessed, "46,XX,del(5)(q13)[10]")
  expect_equal(result$status, "fixed")
})

test_that("fullwidth ISCN punctuation: semicolon, tilde, and equals normalized by preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo(
    "46,XY,t(9\uFF1B22)(q34\uFF1Bq11)[5\uFF5E10]"
  ))
  expect_equal(result$preprocessed, "46,XY,t(9;22)(q34;q11)[5~10]")
  expect_equal(result$status, "fixed")
})

test_that("stray non-ASCII character: no known fix, preprocess_karyo returns NA", {
  result <- suppressMessages(preprocess_karyo(
    "46,XX,der(1)t(1;7)(p11;p11)\u00B5"
  ))
  expect_true(is.na(result$preprocessed))
  expect_equal(result$status, "unfixable")
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

test_that("preprocess_karyo: non-Latin trailing narrative is flagged unfixable (stray_non_ascii), not silently dropped", {
  x <- "46,XY[20] \u30C6\u30B9\u30C8"
  chk <- suppressMessages(check_karyo(x))
  expect_equal(has_issue(x, "trailing_narrative"), 0L)
  expect_equal(chk$stray_non_ascii, 1L)
  expect_equal(chk$unfixable, 1L)
  result <- suppressMessages(preprocess_karyo(x))
  expect_equal(result$preprocessed, NA_character_)
  expect_equal(result$status, "unfixable")
})

test_that("preprocess_karyo: trailing_narrative's fix does not truncate a malformed later clone", {
  result <- suppressMessages(preprocess_karyo("46,XX[10]/46,XX,del(5)(q13)[5"))
  expect_equal(result$preprocessed, NA_character_)
  expect_equal(result$status, "unfixable")
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

test_that("preprocess_karyo: .// and // prefix stripped to donor (fixable)", {
  result <- suppressMessages(preprocess_karyo(c(".//46,XX", "..//46,XX")))
  expect_equal(result$preprocessed, c("46,XX", "46,XX"))
  expect_equal(result$status, c("fixed", "fixed"))
})

test_that("preprocess_karyo: on_chimeric='host' makes zero_host_chimera NA (fixable)", {
  result <- suppressMessages(preprocess_karyo(".//46,XX", on_chimeric = "host"))
  expect_equal(result$preprocessed, NA_character_)
  expect_equal(result$status, "fixed")
})

test_that("preprocess_karyo: on_chimeric='donor' keeps everything after //", {
  result <- suppressMessages(preprocess_karyo(
    "46,XX[15]//47,XY,+8[5]",
    on_chimeric = "donor"
  ))
  expect_equal(result$preprocessed, "47,XY,+8[5]")
  expect_equal(result$status, "fixed")
})

test_that("preprocess_karyo: two or more // are unfixable", {
  result <- suppressMessages(preprocess_karyo("46,XX//47//48"))
  expect_equal(result$preprocessed, NA_character_)
  expect_equal(result$status, "unfixable")
})

test_that("preprocess_karyo: ncSCA token stripped (leading clone)", {
  result <- suppressMessages(preprocess_karyo("ncSCA[4]/46,XY[11]"))
  expect_equal(result$preprocessed, "46,XY[11]")
  expect_equal(result$status, "fixed")
})

test_that("preprocess_karyo: ncSCA token stripped (parenthesized, chimeric host clone)", {
  result <- suppressMessages(preprocess_karyo("46,XX(ncSCA)[1]//46,XY[19]"))
  expect_equal(result$preprocessed, "46,XX[1]")
  expect_equal(result$status, "fixed")
})

test_that("preprocess_karyo: strips leading dot before digit", {
  result <- suppressMessages(preprocess_karyo(c(".46,XX", "..47,XY,+21")))
  expect_equal(result$preprocessed, c("46,XX", "47,XY,+21"))
  expect_equal(result$status, c("fixed", "fixed"))
})

test_that("preprocess_karyo: strips leading dot followed by whitespace before digit", {
  result <- suppressMessages(preprocess_karyo(c(". 47,X,-Y,+1[5]", ".  46,XX")))
  expect_equal(result$preprocessed, c("47,X,-Y,+1[5]", "46,XX"))
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

test_that("preprocess_karyo: empty data frame keeps the detected id column first", {
  df <- data.frame(
    sample_id = character(0),
    karyotype = character(0),
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(
    preprocess_karyo(
      df,
      karyotype_column = "karyotype",
      id_column = "sample_id"
    )
  )
  expect_equal(nrow(result), 0L)
  expect_named(result, c("sample_id", "original", "preprocessed", "status"))
  expect_type(result$sample_id, "character")
})

test_that("preprocess_karyo: rejects input that is neither character, data frame, nor karyo_check", {
  expect_error(
    preprocess_karyo(list("46,XX")),
    "must be a character vector, data frame, or karyo_check object"
  )
  expect_error(preprocess_karyo(1:3), "must be a character vector")
})

test_that("fish_notation: detected when nuc ish suffix present", {
  k <- "[12]/46,XX[3] .nuc ish(PDGFRA x3)[20/200]"
  expect_equal(has_issue(k, "fish_notation"), 1L)
})

test_that("fish_notation: claims trailing FISH clause before trailing_narrative can", {
  k <- "[12]/46,XX[3] .nuc ish(PDGFRA x3)[20/200]"
  expect_equal(has_issue(k, "fish_notation"), 1L)
  expect_equal(has_issue(k, "trailing_narrative"), 0L)

  k2 <- "46,XX,t(9;22)[15] .ish(BCR-ABL)"
  expect_equal(has_issue(k2, "fish_notation"), 1L)
  expect_equal(has_issue(k2, "trailing_narrative"), 0L)
})

test_that("fish_notation: detected when .ish suffix present", {
  k <- "46,XX,t(9;22)[15] .ish(BCR-ABL)"
  expect_equal(has_issue(k, "fish_notation"), 1L)
})

test_that("fish_notation: not detected for clean karyotype", {
  k <- "46,XX,t(9;22)(q34;q11)[15]/46,XX[5]"
  expect_equal(has_issue(k, "fish_notation"), 0L)
})

test_that("preprocess_karyo: strips nuc ish suffix after bracket", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,t(9;22)[15]/46,XX[3] .nuc ish(PDGFRA x3)[20/200]"
    ))$preprocessed,
    "46,XX,t(9;22)[15]/46,XX[3]"
  )
})

test_that("preprocess_karyo: nuc ish block with bracketed FISH counts and interleaved narrative resolves to clean karyotype", {
  k <- paste0(
    "46,XX[20]\nFemale karyotype with no evidence of clonal abnormality\n",
    "nuc ish(IGHx2)[198/200], nuc ish(SNRPN,TP53)x2[199/200]\n",
    "nuc ish(DLEU1,D13S1825)x2[199/200], nuc ish(CDKN2C,CKS1B)x2[200]\n",
    "Normal FISH results for the 14q32 (IGH), 15q, 17p (TP53), 13q loci."
  )
  expect_equal(
    suppressMessages(preprocess_karyo(k))$preprocessed,
    "46,XX[20]"
  )
})

test_that("preprocess_karyo: nuc ish with bare [N] count strips the count (FISH cell count, not metaphase)", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX[20] nuc ish(MYCx2)[200]"
    ))$preprocessed,
    "46,XX[20]"
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
  k <- "46,XX,der(17)t(11;17)(q14;q11).ish der(15)t(15;17)(RARA+,PML+)"
  expect_equal(has_issue(k, "fish_notation"), 1L)
})

test_that("fish_notation: detected for dmin.ish attachment", {
  k <- "46,XY,10~>50dmin.ish del(8)(q24q24)(MYC-),dmin(MYC+)"
  expect_equal(has_issue(k, "fish_notation"), 1L)
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

test_that("fish_notation: detected when .ish appears mid-clone before metaphase count", {
  k <- "46,XX.ish der(10)ins(10;11)(p13;q23q23)(KMT2A+)[20]"
  expect_equal(has_issue(k, "fish_notation"), 1L)
})

test_that("preprocess_karyo: mid-clone .ish stripped, metaphase count preserved", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX.ish der(10)ins(10;11)(p13;q23q23)(KMT2A+)[20]"
    ))$preprocessed,
    "46,XX[20]"
  )
})

test_that("preprocess_karyo: mid-clone .ish stripped per clone in chimeric string", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XY.ish der(10)ins(10;11)(p12;q23q23)(MLL+),der(11)ins(10;11)(p12;q23q23)(MLL+)[2]/47,XY,+8.ish der(10)ins(10;11)(p12;q23q23)(MLL+)[5]/46,XY[13]"
    ))$preprocessed,
    "46,XY[2]/47,XY,+8[5]/46,XY[13]"
  )
})

test_that("preprocess_karyo: mid-clone .ish stripped, conventional aberrations preserved", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XY,der(3)t(3;17)(q29;q21),t(10;11)(p12;q23).ish der(11)t(10;11)(p12;q23)(KMT2A+)[5]/46,XY[16]"
    ))$preprocessed,
    "46,XY,der(3)t(3;17)(q29;q21),t(10;11)(p12;q23)[5]/46,XY[16]"
  )
})

test_that("preprocess_karyo: mid-clone .ish stripped, inv aberration preserved", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,inv(16)(p13q22).ish inv(16)(p13)(MYH11-,CBFB-)(q22)(MYH11+,CBFB+)[20]"
    ))$preprocessed,
    "46,XX,inv(16)(p13q22)[20]"
  )
})

test_that("preprocess_karyo: mid-clone .ish stripped, cp metaphase count preserved", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XY.ish der(10)(PROBE+)[cp10]/47,XY[20]"
    ))$preprocessed,
    "46,XY[cp10]/47,XY[20]"
  )
})

test_that("preprocess_karyo: mid-clone .ish stripped, range metaphase count preserved", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XY.ish der(10)(PROBE+)[1~5]/47,XY[20]"
    ))$preprocessed,
    "46,XY[1~5]/47,XY[20]"
  )
})

test_that("mar_space: detected when space between count and mar", {
  expect_equal(has_issue("47,XY,+1~4 mar[cp15]", "mar_space"), 1L)
})

test_that("mar_space: not detected for well-formed mar token", {
  expect_equal(has_issue("47,XY,+mar[5]", "mar_space"), 0L)
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
  k <- "46,XY,del(5)(q13), .+8[10]/46,XY[5]"
  expect_equal(has_issue(k, "midstring_linewrap"), 1L)
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
