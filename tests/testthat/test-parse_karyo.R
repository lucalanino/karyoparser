test_that("explicit monosomy -7 detected", {
  r <- pk("45,XY,-7")
  expect_equal(r$mono7, 1L)
})

test_that("explicit trisomy +8 detected", {
  r <- pk("47,XY,+8")
  expect_equal(r$tris8, 1L)
})

test_that("-X and +Y detected", {
  r <- pk("46,X,-X,+Y")
  expect_equal(r$monoX, 1L)
  expect_equal(r$trisY, 1L)
})

test_that("del(7) does NOT trigger mono7", {
  r <- pk("46,XY,del(7)(q22)")
  expect_equal(r$mono7, 0L)
})

test_that("multiple monosomies detected", {
  r <- pk("43,XY,-5,-7,-17")
  expect_equal(r$mono5, 1L)
  expect_equal(r$mono7, 1L)
  expect_equal(r$mono17, 1L)
})

test_that("complex karyotype threshold at 3 unique aberrations", {
  r2 <- pk("46,XX,del(5)(q13),del(7)(q22)")
  expect_equal(r2$complex_karyotype, 0L)
  r3 <- pk("46,XX,del(5)(q13),del(7)(q22),+8")
  expect_equal(r3$complex_karyotype, 1L)
})

test_that("idem and sl excluded from complex count", {
  r <- pk("46,XX,del(5)(q13),del(7)(q22)[10]/46,idem,+8[5]")
  expect_equal(r$complex_karyotype, 1L)
})

test_that("repeated aberrations across clones are deduped for complex", {
  # same single aberration in 3 clones -> 1 distinct, not complex
  r1 <- pk(
    "46,XX,t(9;22)(q34;q11)[10]/46,XX,t(9;22)(q34;q11)[5]/46,XX,t(9;22)(q34;q11)[3]"
  )
  expect_equal(r1$complex_karyotype, 0L)
  # two aberrations repeated across clones -> 2 distinct, not complex
  r2 <- pk("46,XX,del(5q),+8[10]/46,XX,del(5q),+8[5]/46,XX,del(5q),+8[3]")
  expect_equal(r2$complex_karyotype, 0L)
  # three distinct aberrations spread one-per-clone -> complex
  r3 <- pk(
    "46,XX,del(5q)[10]/46,XX,t(9;22)(q34;q11)[5]/46,XX,inv(16)(p13q22)[3]"
  )
  expect_equal(r3$complex_karyotype, 1L)
})

test_that("repeated monosomies across clones are deduped for monosomal", {
  # same single monosomy in 2 clones -> 1 distinct, not monosomal
  r1 <- pk("45,XX,-7[10]/45,XX,-7[5]")
  expect_equal(r1$monosomal_karyotype, 0L)
  # one monosomy + one structural, both repeated -> monosomal
  r2 <- pk("46,XX,-7,t(9;22)(q34;q11)[10]/46,XX,-7,t(9;22)(q34;q11)[5]")
  expect_equal(r2$monosomal_karyotype, 1L)
  # two distinct monosomies split one-per-clone -> monosomal
  r3 <- pk("45,XX,-7[10]/45,XX,-5[5]")
  expect_equal(r3$monosomal_karyotype, 1L)
})

test_that("2 autosomal monosomies -> monosomal", {
  r <- pk("44,XY,-5,-7")
  expect_equal(r$monosomal_karyotype, 1L)
})

test_that("1 autosomal mono + 1 structural -> monosomal", {
  r <- pk("45,XY,-7,del(5)(q13)")
  expect_equal(r$monosomal_karyotype, 1L)
})

test_that("autosomal mono + general_translocation -> monosomal", {
  r <- pk("45,XX,-5,t(2;8)(p11;q22)")
  expect_equal(r$general_translocation, 1L)
  expect_equal(r$monosomal_karyotype, 1L)
})

test_that("derivative_chromosome detected", {
  r <- pk("46,XX,der(1)t(2;3)(p11;q22)")
  expect_equal(r$general_derivative, 1L)
})

test_that("+der() detected as derivative_chromosome", {
  r <- pk("47,XX,+der(21)")
  expect_equal(r$general_derivative, 1L)
})

test_that("autosomal mono + +der() -> monosomal", {
  r <- pk("45,XX,-7,+der(21)")
  expect_equal(r$general_derivative, 1L)
  expect_equal(r$monosomal_karyotype, 1L)
})

test_that("derivative_chromosome: autosomal mono + der() -> monosomal", {
  r <- pk("45,XX,-7,der(1)t(2;3)(p11;q22)")
  expect_equal(r$general_derivative, 1L)
  expect_equal(r$monosomal_karyotype, 1L)
})

test_that("derivative_chromosome: co-fires with specific rule for same token", {
  r <- pk("46,XX,der(5)t(5;8)(q11;q11)")
  expect_equal(r$t_5q, 1L)
  expect_equal(r$general_derivative, 1L)
})

test_that("derivative_chromosome: CBF-AML der co-fires but monosomal is overridden", {
  r <- pk("45,XX,-7,der(8)t(8;21)(q22;q22)")
  expect_equal(r$t_8_21_q22_q22, 1L)
  expect_equal(r$general_derivative, 1L)
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("ider() detected as derivative_chromosome", {
  r <- pk("46,XX,ider(1)(q10)t(1;4)(p22;q11)")
  expect_equal(r$general_derivative, 1L)
})

test_that("balanced: bare t() sets balanced, not unbalanced", {
  r <- pk("46,XX,t(9;22)(q34;q11)")
  expect_equal(r$balanced_translocation, 1L)
  expect_equal(r$unbalanced_translocation, 0L)
})

test_that("unbalanced: lone der()t() sets unbalanced, fusion flag still fires", {
  r <- pk("46,XY,der(9)t(9;22)(q34;q11)")
  expect_equal(r$unbalanced_translocation, 1L)
  expect_equal(r$balanced_translocation, 0L)
  expect_equal(r$t_9_22_q34_q11, 1L)
  expect_equal(r$general_derivative, 1L)
})

test_that("balanced: reciprocal der pair is balanced, not unbalanced", {
  r <- pk("46,XX,der(5)t(5;17)(q11;q11),der(17)t(5;17)(q11;q11)")
  expect_equal(r$balanced_translocation, 1L)
  expect_equal(r$unbalanced_translocation, 0L)
})

test_that("balanced: reciprocal pair pairs across partner-order swap", {
  r <- pk("46,XX,der(5)t(5;17)(q11;q11),der(17)t(17;5)(q11;q11)")
  expect_equal(r$balanced_translocation, 1L)
  expect_equal(r$unbalanced_translocation, 0L)
})

test_that("unbalanced: whole-arm der(a;b) is unbalanced", {
  r <- pk("46,XX,der(1;7)(q10;p10)")
  expect_equal(r$unbalanced_translocation, 1L)
  expect_equal(r$balanced_translocation, 0L)
})

test_that("mixed: bare t() plus lone der of same t fires both flags", {
  r <- pk("46,XX,t(8;21)(q22;q22),der(8)t(8;21)(q22;q22)")
  expect_equal(r$balanced_translocation, 1L)
  expect_equal(r$unbalanced_translocation, 1L)
})

test_that("no translocation: both balance flags 0", {
  r <- pk("47,XY,+8")
  expect_equal(r$balanced_translocation, 0L)
  expect_equal(r$unbalanced_translocation, 0L)
})

test_that("unbalanced: lone homologous der is unbalanced, recurrent flag still fires", {
  r3 <- pk("46,XX,der(3)t(3;3)(q21;q26)")
  expect_equal(r3$unbalanced_translocation, 1L)
  expect_equal(r3$balanced_translocation, 0L)
  expect_equal(r3$t_3_3_q21_q26, 1L)
  r16 <- pk("46,XX,der(16)t(16;16)(p13;q22)")
  expect_equal(r16$unbalanced_translocation, 1L)
  expect_equal(r16$balanced_translocation, 0L)
  expect_equal(r16$t_16_16_p13_q22, 1L)
})

test_that("balanced: two homologous der of same signature form a reciprocal pair", {
  r <- pk("46,XX,der(16)t(16;16)(p13;q22),der(16)t(16;16)(p13;q22)")
  expect_equal(r$balanced_translocation, 1L)
  expect_equal(r$unbalanced_translocation, 0L)
})

test_that("partial loss: der(5)t(5;17) implies unbal_partial_loss_5q + unbal_partial_loss_17p", {
  r <- pk("46,XY,der(5)t(5;17)(q11;q11)")
  expect_equal(r$unbal_partial_loss_5q, 1L)
  expect_equal(r$unbal_partial_loss_17p, 1L)
  expect_equal(r$unbal_partial_loss, 1L)
  expect_equal(r$del_5q, 0L)
})

test_that("partial loss: real del(5q) does not set any unbal_partial_loss column", {
  r <- pk("46,XX,del(5q)")
  expect_equal(r$del_5q, 1L)
  expect_equal(r$unbal_partial_loss_5q, 0L)
  expect_equal(r$unbal_partial_loss, 0L)
})

test_that("partial loss: balanced reciprocal pair implies no loss", {
  r <- pk("46,XX,der(5)t(5;17)(q11;q11),der(17)t(5;17)(q11;q11)")
  expect_equal(r$unbal_partial_loss, 0L)
  expect_equal(r$unbal_partial_loss_5q, 0L)
  expect_equal(r$unbal_partial_loss_17p, 0L)
})

test_that("partial loss: all arms exposed, including non-myeloid ones", {
  r <- pk("45,XX,der(9)t(9;22)(q34;q11)")
  expect_equal(r$unbal_partial_loss, 1L)
  expect_equal(r$unbal_partial_loss_9q, 1L)
  expect_equal(r$unbal_partial_loss_22p, 1L)
  expect_equal(r$unbal_partial_loss_5q, 0L)
})

test_that("partial loss: sex-chromosome arms are exposed", {
  r <- pk("46,Y,der(X)t(X;5)(q21;q31)")
  expect_equal(r$unbal_partial_loss_Xq, 1L)
  expect_equal(r$unbal_partial_loss_5p, 1L)
  expect_equal(r$unbal_partial_loss, 1L)
})

test_that("partial loss: der named after a non-partner chromosome derives no loss", {
  r <- pk("46,XX,der(8)t(9;22)(q34;q11)")
  expect_equal(r$unbalanced_translocation, 1L)
  expect_equal(r$unbal_partial_loss, 0L)
})

test_that("partial loss: multi-junction der is unbalanced but derives no loss", {
  r <- pk("46,XX,der(22)t(9;22)(q34;q11)t(11;22)(q23;q11)")
  expect_equal(r$unbalanced_translocation, 1L)
  expect_equal(r$unbal_partial_loss, 0L)
})

test_that("three-way: lone der of a three-way t() is unbalanced, no loss", {
  r <- pk("46,XY,der(9)t(9;22;11)(q34;q11;q23)")
  expect_equal(r$unbalanced_translocation, 1L)
  expect_equal(r$balanced_translocation, 0L)
  expect_equal(r$unbal_partial_loss, 0L)
  expect_equal(r$general_derivative, 1L)
})

test_that("three-way: complete bare t(a;b;c) is balanced", {
  r <- pk("46,XY,t(9;22;11)(q34;q11;q23)")
  expect_equal(r$balanced_translocation, 1L)
  expect_equal(r$unbalanced_translocation, 0L)
})

test_that("three-way: complete der-by-der reciprocal set is balanced", {
  r <- pk(paste0(
    "46,XY,der(9)t(9;22;11)(q34;q11;q23),",
    "der(22)t(9;22;11)(q34;q11;q23),",
    "der(11)t(9;22;11)(q34;q11;q23)"
  ))
  expect_equal(r$balanced_translocation, 1L)
  expect_equal(r$unbalanced_translocation, 0L)
  expect_equal(r$unbal_partial_loss, 0L)
})

test_that("three-way: incomplete set (2 of 3 ders) is unbalanced", {
  r <- pk(paste0(
    "46,XY,der(9)t(9;22;11)(q34;q11;q23),der(22)t(9;22;11)(q34;q11;q23)"
  ))
  expect_equal(r$unbalanced_translocation, 1L)
  expect_equal(r$balanced_translocation, 0L)
})

test_that("partial loss: missing breakpoints yield no derived loss", {
  r <- pk("46,XX,der(5)t(5;17)")
  expect_equal(r$unbalanced_translocation, 1L)
  expect_equal(r$unbal_partial_loss, 0L)
})

test_that("differential: general bare t() is balanced, der()t() is not general tx", {
  bare <- pk("46,XX,t(2;14)(q11;q11)")
  expect_equal(bare$general_translocation, 1L)
  expect_equal(bare$balanced_translocation, 1L)
  der <- pk("46,XX,der(2)t(2;14)(q11;q11)")
  expect_equal(der$general_translocation, 0L)
  expect_equal(der$unbalanced_translocation, 1L)
  expect_equal(der$general_derivative, 1L)
})

test_that("1 autosomal mono alone -> not monosomal", {
  r <- pk("45,XY,-7")
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("sex chromosome monosomy alone doesn't qualify", {
  r <- pk("45,X,-Y")
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("marker does NOT count for monosomal", {
  r <- pk("46,XY,-7,+mar")
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("isochromosome counts as structural for monosomal (no counts_for_monosomal)", {
  # monosomal structural presence now comes from general flags, not a rule flag
  expect_equal(pk("45,XX,-7,i(17)(q10)")$monosomal_karyotype, 1L)
  expect_equal(pk("45,XX,-7,i(7)(q10)")$monosomal_karyotype, 1L)
})

test_that("sex chromosome monosomy + structural does NOT qualify as monosomal", {
  r <- pk("45,X,del(5)(q13)")
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("CBF-AML inv(16) overrides monosomal karyotype to 0", {
  r <- pk("45,XX,-7,inv(16)(p13q22)")
  expect_equal(r$inv_16_p13q22, 1L)
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("CBF-AML t(16;16) overrides monosomal karyotype to 0", {
  r <- pk("45,XX,-7,t(16;16)(p13;q22)")
  expect_equal(r$t_16_16_p13_q22, 1L)
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("idem inherits stemline comma count", {
  r <- pk("46,XX,del(5)(q13),del(7)(q22)[10]/46,idem,+8[5]")
  expect_equal(r$comma_count_aberrations, 3L)
})

test_that("clone without idem uses raw count", {
  r <- pk("46,XX,del(5)(q13)")
  expect_equal(r$comma_count_aberrations, 1L)
})

test_that("idem in clone 2 inheriting 0 aberrations from clean stemline", {
  result <- suppressMessages(check_karyo("46,XX[10]/46,idem[5]"))
  expect_equal(result$invalid_idem, 0L)
  r <- pk("46,XX[10]/46,idem[5]")
  expect_false(is.na(r$chromosome_count))
  expect_equal(r$comma_count_aberrations, 0L)
  expect_equal(r$complex_karyotype, 0L)
})

test_that("on_issues='warn': issue rows get NA output", {
  dirty <- c("46,XX", ".47,XY,+21")
  result <- suppressWarnings(pk(dirty))
  expect_equal(nrow(result), 2L)
  expect_false(is.na(result$chromosome_count[1]))
  expect_true(is.na(result$chromosome_count[2]))
})

test_that("on_issues='warn': emits warning for dirty input", {
  dirty <- c("46,XX", ".47,XY,+21")
  expect_warning(
    suppressMessages(parse_karyo(dirty, on_issues = "warn", verbose = FALSE)),
    regexp = "had issues"
  )
})

test_that("on_issues='stop': throws error for dirty input", {
  dirty <- c("46,XX", ".47,XY,+21")
  expect_error(
    parse_karyo(dirty, on_issues = "stop", verbose = FALSE),
    regexp = "data quality issues"
  )
})

test_that("on_issues='stop': throws error for invalid input", {
  expect_error(
    parse_karyo(c("46,XX", NA), on_issues = "stop", verbose = FALSE),
    regexp = "data quality issues"
  )
})

test_that("on_issues='warn': invalid rows get NA", {
  r <- suppressWarnings(pk(c("46,XX", NA)))
  expect_equal(r$normal_karyotype[1], 1L)
  expect_true(is.na(r$normal_karyotype[2]))
})

test_that("on_issues: clean input produces no guard message", {
  clean <- c("46,XX", "47,XY,+21[10]")
  expect_no_message(
    parse_karyo(clean, on_issues = "warn", verbose = FALSE)
  )
})

test_that("error columns: clean row has both = 0", {
  r <- pk("46,XX")
  expect_equal(r$fixable_error, 0L)
  expect_equal(r$unfixable_error, 0L)
})

test_that("error columns: dirty row under on_issues='warn' has fixable_error=1, unfixable_error=0", {
  r <- suppressWarnings(pk(c("46,XX", ".47,XY,+21")))
  expect_equal(r$fixable_error[1], 0L)
  expect_equal(r$fixable_error[2], 1L)
  expect_equal(r$unfixable_error[2], 0L)
})

test_that("error columns: structural error row has unfixable_error=1, fixable_error=0", {
  r <- suppressWarnings(pk(c("46,XX", NA)))
  expect_equal(r$unfixable_error[1], 0L)
  expect_equal(r$unfixable_error[2], 1L)
  expect_equal(r$fixable_error[2], 0L)
})

test_that("error columns: successfully fixed row has fixable_error=1, unfixable_error=0", {
  r <- parse_karyo(
    c("46,XX", ".47,XY,+21"),
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_equal(r$fixable_error[2], 1L)
  expect_equal(r$unfixable_error[2], 0L)
})

test_that("error columns: present in empty result", {
  r <- pk(character(0))
  expect_true("fixable_error" %in% names(r))
  expect_true("unfixable_error" %in% names(r))
})

test_that("error columns: survive dplyr::filter()", {
  r <- suppressWarnings(pk(c("46,XX", NA, "47,XY,+8")))
  filtered <- dplyr::filter(r, !is.na(normal_karyotype))
  expect_equal(nrow(filtered), 2L)
  expect_true("fixable_error" %in% names(filtered))
  expect_true("unfixable_error" %in% names(filtered))
})

test_that("auto-detect sample_id column", {
  df <- data.frame(sample_id = c("A", "B"), karyotype = c("46,XX", "46,XY"))
  r <- pk(df)
  expect_equal(r$sample_id, c("A", "B"))
  expect_true("sample_id" %in% names(r))
})

test_that("auto-detect patient_id column", {
  df <- data.frame(patient_id = c("P1", "P2"), karyotype = c("46,XX", "46,XY"))
  r <- pk(df)
  expect_equal(r$patient_id, c("P1", "P2"))
})

test_that("explicit id_column parameter", {
  df <- data.frame(my_id = c("X", "Y"), karyotype = c("46,XX", "46,XY"))
  r <- pk(df, id_column = "my_id")
  expect_equal(r$my_id, c("X", "Y"))
  expect_equal(names(r)[1], "my_id")
})

test_that("character vector input has no ID column", {
  r <- pk(c("46,XX", "46,XY"))
  expect_false("sample_id" %in% names(r))
  expect_equal(names(r)[1], "original_karyotype")
})

test_that("parse_karyo: data.frame with all invalid rows and ID column preserves ID in output", {
  df <- data.frame(sample_id = "A", karyotype = "XX,+8")
  r <- suppressWarnings(pk(df))
  expect_equal(r$sample_id, "A")
  expect_true(is.na(r$chromosome_count))
  expect_equal(names(r)[1], "sample_id")
})

test_that("explicit karyotype_column parameter works for non-standard column name", {
  df <- data.frame(my_col = c("46,XX", "46,XY"), stringsAsFactors = FALSE)
  r <- pk(df, karyotype_column = "my_col")
  expect_equal(nrow(r), 2L)
  expect_equal(r$normal_karyotype, c(1L, 1L))
})

test_that("46,XX and 46,XY are normal", {
  r <- pk(c("46,XX", "46,XY"))
  expect_equal(r$normal_karyotype, c(1L, 1L))
})

test_that("abnormal karyotype is not normal", {
  r <- pk("47,XY,+8")
  expect_equal(r$normal_karyotype, 0L)
})

test_that("range karyotype is not normal", {
  r <- pk("45~46,XX[cp20]")
  expect_equal(r$normal_karyotype, 0L)
})

test_that("full pipeline: normal karyotype", {
  r <- pk("46,XX")
  expect_equal(r$chromosome_count, 46L)
  expect_equal(r$normal_karyotype, 1L)
  expect_equal(r$complex_karyotype, 0L)
  expect_equal(r$monosomal_karyotype, 0L)
  expect_equal(r$chromosome_count, 46L)
})

test_that("full pipeline: APL t(15;17)", {
  r <- pk("46,XX,t(15;17)(q24;q21)")
  expect_equal(r$t_15_17_q24_q21, 1L)
  expect_equal(r$chromosome_count, 46L)
  expect_equal(r$normal_karyotype, 0L)
  expect_equal(r$complex_karyotype, 0L)
})

test_that("full pipeline: complex monosomal karyotype", {
  r <- pk("43,XY,-5,-7,del(17)(p13),+8")
  expect_equal(r$mono5, 1L)
  expect_equal(r$mono7, 1L)
  expect_equal(r$tris8, 1L)
  expect_equal(r$del_17p, 1L)
  expect_equal(r$complex_karyotype, 1L)
  expect_equal(r$monosomal_karyotype, 1L)
})

test_that("full pipeline: composite karyotype with metaphases", {
  r <- pk("46,XX,del(5)(q13)[10]/47,XX,del(5)(q13),+8[5]")
  expect_true(!is.na(r$total_metaphases))
  expect_equal(r$total_metaphases, 15L)
  expect_equal(r$del_5q, 1L)
  expect_equal(r$tris8, 1L)
})

test_that("full pipeline: hyperdiploid karyotype", {
  r <- pk("51,XY,+8,+11,+13,+19,+21")
  expect_equal(r$chromosome_count, 51L)
  expect_equal(r$tris8, 1L)
  expect_equal(r$tris21, 1L)
})

test_that("empty input returns empty tibble", {
  r <- pk(character(0))
  expect_equal(nrow(r), 0)
  expect_true("original_karyotype" %in% names(r))
})

test_that("schema is identical across parsed, issue, and empty results", {
  parsed <- pk("46,XX")
  empty <- pk(character(0))
  issue <- suppressWarnings(pk(c("46,XX", "not a karyotype at all")))
  # Guards against drift between the schema's single source of truth and the
  # parsed-row pipeline: all three paths must yield the same columns in the
  # same order.
  expect_identical(names(empty), names(parsed))
  expect_identical(names(issue), names(parsed))
})

test_that(".column_catalog lists exactly the output columns, in order", {
  catalog <- karyoparser:::.column_catalog(myeloid_rules)
  parsed <- pk("46,XX")
  # For a bare character-vector input (no id column) the output columns are
  # exactly the catalog, in catalog order.
  expect_identical(catalog$column, names(parsed))
  expect_setequal(
    unique(catalog$class),
    c(
      "meta",
      "rule",
      "general",
      "aneuploidy",
      "summary",
      "balance",
      "loss",
      "status"
    )
  )
  expect_setequal(catalog$type, c("character", "integer"))
})

test_that("columns = NULL (default) returns the full catalog", {
  catalog <- karyoparser:::.column_catalog(myeloid_rules)
  parsed <- pk("46,XX")
  expect_identical(names(parsed), catalog$column)
})

test_that("columns subsets to the requested classes plus meta/status", {
  parsed <- parse_karyo(
    c("46,XX", "47,XY,+21,t(9;22)(q34;q11.2)"),
    on_issues = "warn",
    verbose = FALSE,
    columns = c("summary", "balance")
  )
  catalog <- karyoparser:::.column_catalog(myeloid_rules)
  expected <- catalog$column[
    catalog$class %in% c("meta", "status", "summary", "balance")
  ]
  expect_identical(names(parsed), expected)
  expect_true("balanced_translocation" %in% names(parsed))
  expect_false("mono21" %in% names(parsed))
  expect_false(any(startsWith(names(parsed), "unbal_partial_loss")))
})

test_that("columns always keeps meta and status classes", {
  parsed <- pk("46,XX", columns = "balance")
  expect_true(all(
    c(
      "original_karyotype",
      "preprocessed_karyotype",
      "chromosome_count",
      "total_metaphases",
      "fixable_error",
      "unfixable_error",
      "chimeric_karyotype",
      "chimeric_clone"
    ) %in%
      names(parsed)
  ))
})

test_that("columns errors on an unknown class", {
  expect_error(
    pk("46,XX", columns = "bogus"),
    "Unknown `columns` value"
  )
})

test_that("columns = 'rule' reflects the flag names of a custom `rules` table", {
  custom_rules <- validate_rules(data.frame(
    flag_name = "my_custom_flag",
    regex = "t\\(9;22\\)\\(q34;q11\\)",
    category = "specific_tx"
  ))
  parsed <- parse_karyo(
    "46,XX,t(9;22)(q34;q11)",
    rules = custom_rules,
    on_issues = "warn",
    verbose = FALSE,
    columns = "rule"
  )
  expect_true("my_custom_flag" %in% names(parsed))
  expect_equal(parsed$my_custom_flag, 1L)
  expect_false("t_15_17_q24_q21" %in% names(parsed))
  expect_false("general_translocation" %in% names(parsed))
})

test_that("single karyotype input works", {
  r <- pk("46,XX")
  expect_equal(nrow(r), 1)
})

test_that("version attribute is set", {
  r <- pk("46,XX")
  expect_equal(
    attr(r, "karyoparser_version"),
    as.character(utils::packageVersion("karyoparser"))
  )
})

test_that("parse_karyo() always returns a tibble", {
  r <- parse_karyo("46,XX", on_issues = "warn", verbose = FALSE)
  expect_true(tibble::is_tibble(r))
})

test_that("comma_count_aberrations for simple karyotype", {
  r <- pk("46,XX")
  expect_equal(r$comma_count_aberrations, 0L)
})

test_that("total_metaphases is NA when no brackets", {
  r <- pk("46,XX")
  expect_true(is.na(r$total_metaphases))
})

test_that("duplicate karyotypes produce same results as unique input", {
  input <- c("46,XX,del(5)(q13)", "47,XY,+8", "46,XX,del(5)(q13)")
  r <- pk(input)
  expect_equal(nrow(r), 3)
  cols <- setdiff(names(r), "original_karyotype")
  expect_equal(r[1, cols], r[3, cols])
})

test_that("dedup preserves row order and all columns", {
  unique_input <- c("46,XX", "47,XY,+8")
  duped_input <- c("46,XX", "47,XY,+8", "46,XX", "47,XY,+8")
  r_unique <- pk(unique_input)
  r_duped <- pk(duped_input)
  expect_equal(nrow(r_duped), 4)
  cols <- setdiff(names(r_unique), c("original_karyotype", "issues"))
  expect_equal(r_duped[1, cols], r_duped[3, cols])
  expect_equal(r_duped[2, cols], r_duped[4, cols])
  expect_equal(r_duped[1, cols], r_unique[1, cols])
  expect_equal(r_duped[2, cols], r_unique[2, cols])
})

test_that("dedup works with data.frame input and ID column", {
  df <- data.frame(
    sample_id = c("A", "B", "C"),
    karyotype = c("46,XX,del(7)(q22)", "46,XY", "46,XX,del(7)(q22)")
  )
  r <- pk(df)
  expect_equal(nrow(r), 3)
  expect_equal(r$sample_id, c("A", "B", "C"))
  expect_equal(r$del_7q[1], r$del_7q[3])
})

test_that("dedup handles mix of valid and invalid karyotypes", {
  input <- c("46,XX", NA, "46,XX", "47,XY,+8", NA)
  r <- suppressWarnings(pk(input))
  expect_equal(nrow(r), 5)
  expect_equal(r$normal_karyotype[1], 1L)
  expect_equal(r$normal_karyotype[3], 1L)
  expect_true(is.na(r$normal_karyotype[2]))
  expect_true(is.na(r$normal_karyotype[5]))
})

test_that("on_issues='preprocess': dirty row is cleaned and parsed (not NA)", {
  r <- parse_karyo(
    c("46,XX", ".47,XY,+21"),
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_equal(nrow(r), 2L)
  expect_equal(r$tris21[2], 1L)
  expect_false(is.na(r$chromosome_count[2]))
})

test_that("on_issues='preprocess': structural errors still produce NA", {
  r <- parse_karyo(
    c("46,XX", "not_a_karyotype"),
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_true(is.na(r$chromosome_count[2]))
})

test_that("on_issues='preprocess': verbose message reports fixed count", {
  suppressMessages(expect_message(
    parse_karyo(
      c("46,XX", ".47,XY,+21"),
      on_issues = "preprocess",
      verbose = TRUE
    ),
    "Fixed"
  ))
})

test_that("on_issues='preprocess': successfully fixed row has fixable_error=1, unfixable_error=0", {
  r <- parse_karyo(
    c("46,XX", ".47,XY,+21"),
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_equal(r$fixable_error[2], 1L)
  expect_equal(r$unfixable_error[2], 0L)
})

test_that("on_issues='preprocess': multiple dirty rows all cleaned", {
  r <- parse_karyo(
    c(".46,XX", ".47,XY,+21"),
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_equal(r$normal_karyotype[1], 1L)
  expect_equal(r$tris21[2], 1L)
})

test_that("on_issues='warn': dirty row returns NA", {
  r <- suppressWarnings(parse_karyo(
    c("46,XX", ".47,XY,+21"),
    on_issues = "warn",
    verbose = FALSE
  ))
  expect_true(is.na(r$chromosome_count[2]))
})

test_that("on_issues='warn': clean row in same batch unaffected", {
  r <- suppressWarnings(parse_karyo(
    c("46,XX", ".47,XY,+21"),
    on_issues = "warn",
    verbose = FALSE
  ))
  expect_equal(r$normal_karyotype[1], 1L)
})

test_that(".dirty_patterns contains all expected keys", {
  expect_setequal(
    names(karyoparser:::.dirty_patterns),
    c(
      "unicode_notation",
      "embedded_newline",
      "html_entities",
      "leading_dot",
      "count_sex_separator",
      "fish_notation",
      "trailing_narrative",
      "midstring_linewrap",
      "missing_sex_comma",
      "mar_space",
      "zero_host_chimera"
    )
  )
})

test_that("preprocess_karyo() output is consistent with .dirty_patterns fix rules", {
  result_lt <- suppressMessages(preprocess_karyo("&lt;46&gt;,XX"))
  expect_equal(result_lt$status, "unfixable")
  expect_equal(result_lt$preprocessed, NA_character_)
  expect_equal(
    suppressMessages(preprocess_karyo(".46,XY,+21[10] .Note"))$preprocessed,
    "46,XY,+21[10]"
  )
  result_zhc <- suppressMessages(preprocess_karyo(".//46,XX"))
  expect_equal(result_zhc$status, "fixed")
  expect_equal(result_zhc$preprocessed, "46,XX")
})

test_that("missing_sex_comma: detected when sex complement followed by space+aberration", {
  result <- suppressMessages(check_karyo("46,XX der(15;17)(q10;q10)[10]"))
  expect_equal(result$missing_sex_comma, 1L)
})

test_that("missing_sex_comma: not detected for well-formed karyotype", {
  result <- suppressMessages(check_karyo("46,XX,der(15;17)(q10;q10)[10]"))
  expect_equal(result$missing_sex_comma, 0L)
})

test_that("missing_sex_comma: preprocess inserts missing comma", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX der(15;17)(q10;q10),+17[10]"
    ))$preprocessed,
    "46,XX,der(15;17)(q10;q10),+17[10]"
  )
})

test_that("missing_sex_comma: full pipeline via on_issues=fix parses correctly", {
  k <- ".46,XX der(15;16)(q10;q10),+16,der(16)t(1;16)(q12;q11.2)[11]/46,XX[4] .Abnormal clone"
  res <- parse_karyo(k, on_issues = "preprocess", verbose = FALSE)
  expect_equal(res$normal_karyotype, 0L)
  expect_equal(res$chromosome_count, 46L)
})

test_that("missing_sex_comma: detected and fixed for XXYY and XXXY", {
  expect_equal(
    suppressMessages(check_karyo(
      "48,XXYY der(5;17)(q10;q10)[10]"
    ))$missing_sex_comma,
    1L
  )
  expect_equal(
    suppressMessages(check_karyo(
      "48,XXXY der(5;17)(q10;q10)[10]"
    ))$missing_sex_comma,
    1L
  )
  expect_equal(
    suppressMessages(preprocess_karyo(
      "48,XXYY der(5;17)(q10;q10)[10]"
    ))$preprocessed,
    "48,XXYY,der(5;17)(q10;q10)[10]"
  )
  expect_equal(
    suppressMessages(preprocess_karyo(
      "48,XXXY der(5;17)(q10;q10)[10]"
    ))$preprocessed,
    "48,XXXY,der(5;17)(q10;q10)[10]"
  )
})

test_that("count_sex_separator: missing comma between count and sex is repaired", {
  result <- suppressMessages(check_karyo(
    "46XY,der(7)t(7;11)(q11.2;q13)[2]/46,XY[9]"
  ))
  expect_equal(result$count_sex_separator, 1L)
  expect_equal(result$fixable, 1L)
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46XY,der(7)t(7;11)(q11.2;q13)[2]/46,XY[9]"
    ))$preprocessed,
    "46,XY,der(7)t(7;11)(q11.2;q13)[2]/46,XY[9]"
  )
})

test_that("count_sex_separator: dot between count and sex is replaced with comma", {
  result <- suppressMessages(check_karyo("45.XY,-7[7]/46,XY[8]"))
  expect_equal(result$count_sex_separator, 1L)
  r <- parse_karyo(
    "45.XY,-7[7]/46,XY[8]",
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_equal(r$preprocessed_karyotype, "45,XY,-7[7]/46,XY[8]")
  expect_equal(r$mono7, 1L)
})

test_that("count_sex_separator: well-formed clones are not touched", {
  for (k in c("46,XX", "47,XY,+21[10]", "46,XX[15]//46,XY[5]")) {
    expect_equal(
      suppressMessages(check_karyo(k))$count_sex_separator,
      0L,
      info = k
    )
  }
})

test_that("preprocess_karyo: anchor strips trailing content after last bracket", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,add(3)(p13),-22,+mar[4] ./46,XX[11] .Abnormal detected"
    ))$preprocessed,
    "46,XX,add(3)(p13),-22,+mar[4]/46,XX[11]"
  )
})

test_that("preprocess_karyo: anchor handles cp bracket", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,del(5)(q13)[cp10] .Some note"
    ))$preprocessed,
    "46,XX,del(5)(q13)[cp10]"
  )
})

test_that("preprocess_karyo: anchor handles range bracket", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46~48,XX[5~20] .Some note"
    ))$preprocessed,
    "46~48,XX[5~20]"
  )
})

test_that("preprocess_karyo: anchor on no-bracket string leaves unchanged", {
  expect_equal(
    suppressMessages(preprocess_karyo("46,XX"))$preprocessed,
    "46,XX"
  )
})

test_that("preprocess_karyo: trailing content with its own bracket — rule 3 cleans up remainder", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX[20] .Note with [5] cells"
    ))$preprocessed,
    "46,XX[20]"
  )
})

test_that("preprocess_karyo: multi-clone clean string unchanged by anchor", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,t(9;22)[15]/46,XX[5]"
    ))$preprocessed,
    "46,XX,t(9;22)[15]/46,XX[5]"
  )
})

test_that("preprocess_karyo: ] ./ separator collapsed to ]/", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,+mar[cp8] ./46,XX[7]"
    ))$preprocessed,
    "46,XX,+mar[cp8]/46,XX[7]"
  )
})

test_that("preprocess_karyo: ] ./ mid-string and trailing narrative both fixed", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,+mar[4] ./46,XX[11] .Abnormal clone detected"
    ))$preprocessed,
    "46,XX,+mar[4]/46,XX[11]"
  )
})

test_that("preprocess_karyo: mid-string ', .der()' collapsed", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,add(12)(p13),-13, .der(13;14)(q10;q10),-14[10]/46,XX[5]"
    ))$preprocessed,
    "46,XX,add(12)(p13),-13,der(13;14)(q10;q10),-14[10]/46,XX[5]"
  )
})

test_that("preprocess_karyo: mid-string ', .del()' collapsed", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,+8, .del(5)(q13q33)[10]/46,XX[5]"
    ))$preprocessed,
    "46,XX,+8,del(5)(q13q33)[10]/46,XX[5]"
  )
})

test_that("preprocess_karyo: uppercase after ', .' not collapsed by midstring_linewrap", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX[20] .Female karyotype"
    ))$preprocessed,
    "46,XX[20]"
  )
})

test_that("check_karyo detects chimeric_separator as fixable", {
  result <- suppressMessages(check_karyo("46,XX[15]//46,XY[5]"))
  expect_equal(result$chimeric_separator, 1L)
  expect_equal(result$fixable, 1L)
})

test_that("check_karyo: clean karyotype does not trigger chimeric_separator", {
  result <- suppressMessages(check_karyo("46,XX[20]"))
  expect_equal(result$chimeric_separator, 0L)
})

test_that("check_karyo detects updated_iscn as unfixable", {
  result <- suppressMessages(check_karyo(
    "46,XX,add(9)[3]/46,XY[12] Updated ISCN 45,XY[15]"
  ))
  expect_equal(result$updated_iscn, 1L)
  expect_equal(result$unfixable, 1L)
})

test_that("check_karyo updated_iscn is case-insensitive", {
  result <- suppressMessages(check_karyo("46,XX updated iscn new version"))
  expect_equal(result$updated_iscn, 1L)
})

test_that("updated_iscn rows become NA in parse_karyo", {
  r <- suppressWarnings(pk("46,XX Updated ISCN new"))
  expect_true(is.na(r$chromosome_count))
})

test_that("updated_iscn rows stay NA even with on_issues='preprocess'", {
  r <- parse_karyo(
    "46,XX Updated ISCN new",
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_true(is.na(r$chromosome_count))
})

test_that("on_issues='warn': chimeric row is still clone-selected and parsed", {
  r <- suppressMessages(suppressWarnings(
    parse_karyo("46,XX[15]//46,XY[5]", on_issues = "warn", verbose = FALSE)
  ))
  expect_false(is.na(r$chromosome_count))
  expect_equal(r$chimeric_clone, "host")
})

test_that("on_issues='warn': chimeric row has fixable_error=1", {
  r <- suppressMessages(suppressWarnings(
    parse_karyo("46,XX[15]//46,XY[5]", on_issues = "warn", verbose = FALSE)
  ))
  expect_equal(r$fixable_error, 1L)
  expect_equal(r$unfixable_error, 0L)
})

test_that("on_issues='preprocess': chimeric row truncated and parsed", {
  r <- suppressMessages(suppressWarnings(parse_karyo(
    "46,XX,t(9;22)(q34;q11)[15]//46,XX[5]",
    on_issues = "preprocess",
    verbose = FALSE
  )))
  expect_false(is.na(r$chromosome_count))
  expect_equal(r$t_9_22_q34_q11, 1L)
})

test_that("on_issues='preprocess': clean row in same batch unaffected by chimeric truncation", {
  r <- suppressMessages(suppressWarnings(parse_karyo(
    c("46,XX[20]", "46,XX,t(9;22)(q34;q11)[15]//46,XX[5]"),
    on_issues = "preprocess",
    verbose = FALSE
  )))
  expect_equal(r$normal_karyotype[1], 1L)
  expect_false(is.na(r$chromosome_count[2]))
})

test_that("on_issues='preprocess': verbose message reports fixed count for chimeric row", {
  suppressMessages(suppressWarnings(expect_message(
    parse_karyo(
      "46,XX,t(9;22)(q34;q11)[15]//46,XX[5]",
      on_issues = "preprocess",
      verbose = TRUE
    ),
    "Fixed"
  )))
})

test_that("on_issues='preprocess': chimeric row with residual structural issue after truncation becomes NA", {
  r <- suppressMessages(parse_karyo(
    "not_valid//46,XX[5]",
    on_issues = "preprocess",
    verbose = FALSE
  ))
  expect_true(is.na(r$chromosome_count))
})

test_that("on_issues='stop': clean chimeric row does not error, is parsed", {
  r <- suppressMessages(
    parse_karyo("46,XX[15]//46,XY[5]", on_issues = "stop", verbose = FALSE)
  )
  expect_false(is.na(r$chromosome_count))
  expect_equal(r$chimeric_clone, "host")
})

test_that("on_issues='stop': non-chimeric issue still errors despite chimeric row", {
  expect_error(
    parse_karyo(
      c("46,XX[15]//46,XY[5]", "not_a_karyotype"),
      on_issues = "stop",
      verbose = FALSE
    ),
    "stopped"
  )
})

test_that("check_karyo detects zero_host_chimera as fixable", {
  result <- suppressMessages(check_karyo(".//46,XX[10]"))
  expect_equal(result$zero_host_chimera, 1L)
  expect_equal(result$unfixable, 0L)
  expect_equal(result$fixable, 1L)
})

test_that("check_karyo detects multiple_chimeric_separator as unfixable", {
  result <- suppressMessages(check_karyo("46,XX//47//48"))
  expect_equal(result$multiple_chimeric_separator, 1L)
  expect_equal(result$unfixable, 1L)
  expect_equal(result$fixable, 0L)
})

test_that("check_karyo: multiple leading dots also detected as zero_host_chimera", {
  result <- suppressMessages(check_karyo("..//46,XY[5]"))
  expect_equal(result$zero_host_chimera, 1L)
})

test_that("check_karyo: normal karyotype does not trigger zero_host_chimera", {
  result <- suppressMessages(check_karyo("46,XX[20]"))
  expect_equal(result$zero_host_chimera, 0L)
})

test_that("parse_karyo: zero_host_chimera parsed to donor under default", {
  r <- suppressMessages(suppressWarnings(pk(".//46,XX[10]")))
  expect_false(is.na(r$chromosome_count))
  expect_equal(r$chimeric_clone, "donor")
})

test_that("parse_karyo: zero_host_chimera has fixable_error=1, unfixable_error=0", {
  r <- suppressMessages(suppressWarnings(pk(".//46,XX[10]")))
  expect_equal(r$unfixable_error, 0L)
  expect_equal(r$fixable_error, 1L)
})

test_that("parse_karyo: zero_host_chimera NA but fixable under on_chimeric='host'", {
  r <- suppressMessages(
    parse_karyo(
      ".//46,XX[10]",
      on_issues = "preprocess",
      on_chimeric = "host",
      verbose = FALSE
    )
  )
  expect_true(is.na(r$chromosome_count))
  expect_equal(r$unfixable_error, 0L)
  expect_equal(r$fixable_error, 1L)
  expect_equal(r$chimeric_karyotype, 1L)
})

test_that("parse_karyo: zero_host_chimera with other dirty patterns is parsed", {
  r <- suppressMessages(parse_karyo(
    ".//46,XX[10] .Female karyotype",
    on_issues = "preprocess",
    verbose = FALSE
  ))
  expect_false(is.na(r$chromosome_count))
  expect_equal(r$chimeric_clone, "donor")
})

test_that("parse_karyo: zero_host + trailing narrative is fully fixable", {
  r <- suppressMessages(suppressWarnings(pk(".//46,XX[10] .Female karyotype")))
  expect_equal(r$fixable_error, 1L)
  expect_equal(r$unfixable_error, 0L)
})

test_that("check_karyo: verbose=TRUE prints checking count and summary", {
  suppressMessages(expect_message(
    check_karyo(c("46,XX", ".47,XY,+21", NA), verbose = TRUE),
    "Checking 3 karyotype"
  ))
  suppressMessages(expect_message(
    check_karyo(c("46,XX", ".47,XY,+21", NA), verbose = TRUE),
    "Fixable"
  ))
})

test_that("check_karyo: verbose=TRUE prints 'All clean.' when no issues", {
  suppressMessages(expect_message(
    check_karyo("46,XX[20]", verbose = TRUE),
    "All clean"
  ))
})

test_that("check_karyo: verbose=TRUE adds per-type breakdown", {
  suppressMessages(expect_message(
    check_karyo(c("46,XX", ".47,XY,+21", NA), verbose = TRUE),
    "breakdown"
  ))
})

test_that("parse_karyo: clean row in same batch unaffected by zero_host_chimera row", {
  r <- suppressMessages(suppressWarnings(pk(c("46,XX[20]", ".//46,XY[10]"))))
  expect_equal(r$normal_karyotype[1], 1L)
  expect_false(is.na(r$chromosome_count[2]))
  expect_equal(r$chimeric_clone, c(NA_character_, "donor"))
})

test_that("midstring_linewrap: dot variant still detected", {
  result <- suppressMessages(check_karyo(
    "46,XY,del(5)(q13), .t(9;22)(q34;q11.2)[10]/46,XY[5]"
  ))
  expect_equal(result$midstring_linewrap, 1L)
})

test_that("preprocess_karyo: collapses ', .t()' (dot present) mid-string artifact", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XY,del(5)(q13), .t(9;22)(q34;q11.2)[10]/46,XY[5]"
    ))$preprocessed,
    "46,XY,del(5)(q13),t(9;22)(q34;q11.2)[10]/46,XY[5]"
  )
})

test_that("normalize_iscn: bare ', t()' space (no dot) cleaned via preprocess_karyo", {
  result <- suppressMessages(
    preprocess_karyo("46,XY,del(5)(q13), t(9;22)(q34;q11.2)[10]/46,XY[5]")
  )
  expect_equal(
    result$preprocessed,
    "46,XY,del(5)(q13),t(9;22)(q34;q11.2)[10]/46,XY[5]"
  )
})

test_that("trailing_narrative: detected when narrative follows last paren (no bracket)", {
  result <- suppressMessages(check_karyo(
    "46,XX,t(9;22)(q34;q11.2) Abnormal female karyotype"
  ))
  expect_equal(result$trailing_narrative, 1L)
})

test_that("preprocess_karyo: strips narrative after last ')' when no bracket present", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,t(9;22)(q34;q11.2) Abnormal female karyotype"
    ))$preprocessed,
    "46,XX,t(9;22)(q34;q11.2)"
  )
})

test_that("preprocess_karyo: trailing narrative rule 3 does not fire when bracket present", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,t(9;22)(q34;q11.2)[20] Abnormal female karyotype"
    ))$preprocessed,
    "46,XX,t(9;22)(q34;q11.2)[20]"
  )
})

test_that("trailing_narrative: detected when Capital follows bracket without dot", {
  result <- suppressMessages(check_karyo("46,XX[20] Female karyotype"))
  expect_equal(result$trailing_narrative, 1L)
  expect_equal(result$fixable, 1L)
})

test_that("preprocess_karyo: strips narrative after bracket when no dot present", {
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX[20] Female karyotype"
    ))$preprocessed,
    "46,XX[20]"
  )
})

test_that("trailing_narrative: stripped after count+sex with no bracket or paren", {
  for (k in c(
    "46,XY Normal male karyotype",
    "46,XX  Normal female karyotype",
    "45~46,XY Some narrative"
  )) {
    chk <- suppressMessages(check_karyo(k))
    expect_equal(chk$trailing_narrative, 1L, info = k)
    expect_equal(chk$no_sex_complement, 0L, info = k)
    expect_equal(chk$fixable, 1L, info = k)
  }
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XY Normal male karyotype"
    ))$preprocessed,
    "46,XY"
  )
  # A bare count+sex with no narrative stays clean (rule does not over-fire).
  expect_equal(suppressMessages(check_karyo("46,XY"))$fixable, 0L)
})

test_that("trailing_narrative rule 5 does not strip the Updated ISCN marker", {
  chk <- suppressMessages(check_karyo(
    "46,XX,add(9)[3]/46,XY[12] Updated ISCN 45,XY[15]"
  ))
  expect_equal(chk$updated_iscn, 1L)
  expect_equal(chk$unfixable, 1L)
})

test_that("normalize_iscn: collapses space before opening bracket", {
  result <- suppressMessages(pk("46,XX [20]"))
  expect_equal(result$normal_karyotype, 1L)
})

test_that("preprocess_karyo: space before '[' cleaned via normalize_iscn in preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo("46,XX [20]"))
  expect_equal(result$preprocessed, "46,XX[20]")
})

test_that("parse_karyo: chimeric_karyotype column exists in output", {
  result <- suppressMessages(pk("46,XX"))
  expect_true("chimeric_karyotype" %in% names(result))
})

test_that("parse_karyo: clean row has chimeric_karyotype=0", {
  result <- suppressMessages(pk("46,XX"))
  expect_equal(result$chimeric_karyotype, 0L)
})

test_that("parse_karyo on_issues='preprocess': chimeric row has chimeric_karyotype=1", {
  result <- suppressWarnings(suppressMessages(
    parse_karyo(
      "46,XX,t(9;22)(q34;q11.2)[3]//46,XY[12]",
      on_issues = "preprocess",
      verbose = FALSE
    )
  ))
  expect_equal(result$chimeric_karyotype, 1L)
})

test_that("parse_karyo on_issues='warn': chimeric row has chimeric_karyotype=1", {
  result <- suppressWarnings(suppressMessages(
    parse_karyo(
      "46,XX,t(9;22)(q34;q11.2)[3]//46,XY[12]",
      on_issues = "warn",
      verbose = FALSE
    )
  ))
  expect_equal(result$chimeric_karyotype, 1L)
})

test_that("parse_karyo: zero_host_chimera row has chimeric_karyotype=1", {
  result <- suppressMessages(
    parse_karyo(".//46,XX[20]", on_issues = "preprocess", verbose = FALSE)
  )
  expect_equal(result$chimeric_karyotype, 1L)
})

test_that("parse_karyo: zero_host_chimera has fixable_error=1 and chimeric_karyotype=1", {
  result <- suppressMessages(
    parse_karyo(".//46,XX[20]", on_issues = "preprocess", verbose = FALSE)
  )
  expect_equal(result$fixable_error, 1L)
  expect_equal(result$chimeric_karyotype, 1L)
})

test_that("parse_karyo: multiple_chimeric_separator has unfixable_error=1, chimeric_karyotype=1", {
  result <- suppressMessages(
    parse_karyo("46,XX//47//48", on_issues = "preprocess", verbose = FALSE)
  )
  expect_true(is.na(result$chromosome_count))
  expect_equal(result$unfixable_error, 1L)
  expect_equal(result$chimeric_karyotype, 1L)
  expect_true(is.na(result$chimeric_clone))
})

test_that("parse_karyo: non-chimeric row in chimeric batch has chimeric_karyotype=0", {
  result <- suppressWarnings(suppressMessages(
    parse_karyo(
      c("46,XX,t(9;22)(q34;q11.2)[3]//46,XY[12]", "46,XY"),
      on_issues = "preprocess",
      verbose = FALSE
    )
  ))
  expect_equal(result$chimeric_karyotype, c(1L, 0L))
})

test_that("parse_karyo: chimeric_clone column exists and is NA for clean rows", {
  result <- suppressMessages(pk("46,XX"))
  expect_true("chimeric_clone" %in% names(result))
  expect_true(is.na(result$chimeric_clone))
})

test_that("parse_karyo: on_chimeric='donor' keeps donor population after //", {
  result <- suppressMessages(suppressWarnings(parse_karyo(
    "46,XX[15]//47,XY,+8[5]",
    on_issues = "preprocess",
    on_chimeric = "donor"
  )))
  expect_equal(result$preprocessed_karyotype, "47,XY,+8[5]")
  expect_equal(result$chimeric_clone, "donor")
  expect_equal(result$tris8, 1L)
})

test_that("parse_karyo: on_chimeric='host' keeps host, default keeps host too", {
  res_host <- suppressMessages(suppressWarnings(parse_karyo(
    "47,XY,+8[15]//46,XX[5]",
    on_issues = "preprocess",
    on_chimeric = "host"
  )))
  expect_equal(res_host$chimeric_clone, "host")
  expect_equal(res_host$tris8, 1L)
})

test_that("parse_karyo: on_chimeric inherited from karyo_preprocessed", {
  pp <- suppressMessages(preprocess_karyo(
    "46,XX[15]//47,XY,+8[5]",
    on_chimeric = "donor"
  ))
  result <- suppressMessages(suppressWarnings(parse_karyo(pp)))
  expect_equal(result$chimeric_clone, "donor")
  expect_equal(result$tris8, 1L)
})

test_that("parse_karyo: conflicting explicit on_chimeric on preprocessed input errors", {
  pp <- suppressMessages(preprocess_karyo(
    "46,XX[15]//47,XY,+8[5]",
    on_chimeric = "donor"
  ))
  expect_error(
    suppressMessages(parse_karyo(pp, on_chimeric = "host")),
    "conflicts with the karyo_preprocessed input"
  )
})

test_that("preprocessed_karyotype: normalize_iscn applied even in on_issues='warn'", {
  r <- parse_karyo("46 , XX", on_issues = "warn", verbose = FALSE)
  expect_equal(r$preprocessed_karyotype, "46,XX")
})

test_that("preprocessed_karyotype: consistent across all three on_issues modes for clean input", {
  x <- "46,XX,t(9;22)(q34;q11.2)[20]"
  r_fix <- parse_karyo(x, on_issues = "preprocess", verbose = FALSE)
  r_warn <- parse_karyo(x, on_issues = "warn", verbose = FALSE)
  expect_equal(r_fix$preprocessed_karyotype, x)
  expect_equal(r_warn$preprocessed_karyotype, x)
})

test_that("fixable_error=1 for a fixed row in on_issues='fix' mode", {
  r <- parse_karyo(
    ".46,XX,t(9;22)(q34;q11.2)[20]",
    on_issues = "preprocess",
    verbose = FALSE
  )
  expect_equal(r$fixable_error, 1L)
  expect_equal(r$unfixable_error, 0L)
})

test_that("fixable_error consistent: same fixable row gives fixable_error=1 in both modes", {
  x <- ".46,XX,t(9;22)(q34;q11.2)[20]"
  r_fix <- suppressWarnings(parse_karyo(
    x,
    on_issues = "preprocess",
    verbose = FALSE
  ))
  r_warn <- suppressWarnings(parse_karyo(
    x,
    on_issues = "warn",
    verbose = FALSE
  ))
  expect_equal(r_fix$fixable_error, 1L)
  expect_equal(r_warn$fixable_error, 1L)
  expect_false(is.na(r_fix$chromosome_count))
  expect_true(is.na(r_warn$chromosome_count))
})

test_that("parse_karyo on_issues='warn': missing_sex_comma does not false-positive no_sex_complement", {
  r <- suppressWarnings(parse_karyo(
    "46,XX der(7)t(7;12)(q36;q24)[10]",
    on_issues = "warn",
    verbose = FALSE
  ))
  expect_equal(r$fixable_error, 1L)
  expect_equal(r$unfixable_error, 0L)
})

test_that("check_karyo: missing_sex_comma does not also fire no_sex_complement", {
  result <- suppressMessages(
    check_karyo("46,XX der(7)t(7;12)(q36;q24)[10]")
  )
  expect_equal(result$missing_sex_comma, 1L)
  expect_equal(result$no_sex_complement, 0L)
})

test_that("preprocess -> parse workflow: zero_host_chimera parsed to donor", {
  raw <- c("46,XX[20]", ".//46,XY[10]", "46,XX,t(9;22)(q34;q11)[15]")
  proc <- suppressMessages(preprocess_karyo(raw))
  expect_equal(
    proc$preprocessed,
    c("46,XX[20]", "46,XY[10]", "46,XX,t(9;22)(q34;q11)[15]")
  )
  expect_equal(proc$status, c("clean", "fixed", "clean"))
  parsed <- suppressMessages(
    parse_karyo(proc, on_issues = "warn")
  )
  expect_false(is.na(parsed$chromosome_count[1]))
  expect_false(is.na(parsed$chromosome_count[2]))
  expect_false(is.na(parsed$chromosome_count[3]))
  expect_equal(parsed$chimeric_clone[2], "donor")
  expect_equal(parsed$fixable_error[2], 1L)
})

test_that("preprocess -> parse workflow: fixable rows are parsed, unfixable are NA", {
  raw <- c(".46,XX[20]", "not_a_karyotype", "46,XY,+8[10]")
  proc <- suppressMessages(preprocess_karyo(raw))
  expect_equal(proc$status, c("fixed", "unfixable", "clean"))
  expect_equal(proc$preprocessed[2], NA_character_)
  parsed <- suppressMessages(
    parse_karyo(proc, on_issues = "warn")
  )
  expect_false(is.na(parsed$chromosome_count[1]))
  expect_true(is.na(parsed$chromosome_count[2]))
  expect_equal(parsed$unfixable_error[2], 1L)
  expect_false(is.na(parsed$chromosome_count[3]))
})

test_that("check_karyo: accepts data frame, auto-detects karyotype column", {
  df <- data.frame(
    karyotype = c("46,XX", ".46,XY[20]"),
    stringsAsFactors = FALSE
  )
  result <- check_karyo(df)
  expect_s3_class(result, "karyo_check")
  expect_equal(nrow(result), 2L)
  expect_equal(result$karyotype, c("46,XX", ".46,XY[20]"))
  expect_equal(result$leading_dot, c(0L, 1L))
})

test_that("check_karyo: data frame with id column propagates id as first column", {
  df <- data.frame(
    sample_id = c("S1", "S2"),
    karyotype = c("46,XX", "46,XY"),
    stringsAsFactors = FALSE
  )
  result <- check_karyo(df)
  expect_equal(names(result)[1], "sample_id")
  expect_equal(result$sample_id, c("S1", "S2"))
})

test_that("check_karyo: iscn column name auto-detected", {
  df <- data.frame(iscn = c("46,XX", "47,XY,+21"), stringsAsFactors = FALSE)
  result <- check_karyo(df)
  expect_s3_class(result, "karyo_check")
  expect_equal(nrow(result), 2L)
})

test_that("check_karyo: explicit karyotype_column", {
  df <- data.frame(my_col = c("46,XX", "46,XY"), stringsAsFactors = FALSE)
  result <- check_karyo(df, karyotype_column = "my_col")
  expect_equal(nrow(result), 2L)
})

test_that("check_karyo: explicit id_column", {
  df <- data.frame(
    my_id = c("A", "B"),
    karyotype = c("46,XX", "46,XY"),
    stringsAsFactors = FALSE
  )
  result <- check_karyo(df, id_column = "my_id")
  expect_equal(names(result)[1], "my_id")
  expect_equal(result$my_id, c("A", "B"))
})

test_that("preprocess_karyo: accepts data frame, auto-detects karyotype column", {
  df <- data.frame(karyotype = c("46,XX", ".46,XY"), stringsAsFactors = FALSE)
  result <- suppressMessages(preprocess_karyo(df, verbose = TRUE))
  expect_s3_class(result, "karyo_preprocessed")
  expect_equal(nrow(result), 2L)
  expect_equal(result$status, c("clean", "fixed"))
})

test_that("preprocess_karyo: data frame with id column propagates id as first column", {
  df <- data.frame(
    sample_id = c("S1", "S2"),
    karyotype = c("46,XX", "46,XY"),
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(preprocess_karyo(df))
  expect_equal(names(result)[1], "sample_id")
  expect_equal(result$sample_id, c("S1", "S2"))
})

test_that("preprocess_karyo: explicit karyotype_column for non-candidate name", {
  df <- data.frame(my_iscn = c("46,XX", "46,XY"), stringsAsFactors = FALSE)
  result <- suppressMessages(preprocess_karyo(df, karyotype_column = "my_iscn"))
  expect_equal(nrow(result), 2L)
  expect_equal(result$status, c("clean", "clean"))
})

test_that("full pipeline check -> preprocess -> parse: no extra args needed", {
  df <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    karyotype = c("46,XX[20]", ".46,XY,+8[10]", ".//46,XX[10]"),
    stringsAsFactors = FALSE
  )
  ck <- check_karyo(df)
  pp <- suppressMessages(preprocess_karyo(ck))
  result <- suppressWarnings(suppressMessages(parse_karyo(pp)))
  expect_true("sample_id" %in% names(result))
  expect_equal(result$sample_id, c("S1", "S2", "S3"))
  expect_false(is.na(result$chromosome_count[1]))
  expect_false(is.na(result$chromosome_count[2]))
  expect_equal(result$fixable_error[2], 1L)
  expect_false(is.na(result$chromosome_count[3]))
  expect_equal(result$chimeric_clone[3], "donor")
  expect_equal(result$fixable_error[3], 1L)
})

test_that("karyo_preprocessed with unfixable rows: default on_issues='stop' warns, not errors", {
  pp <- suppressMessages(preprocess_karyo(c("46,XX", NA)))
  expect_warning(
    result <- parse_karyo(pp),
    regexp = "unfixable"
  )
  expect_false(is.na(result$chromosome_count[1]))
  expect_true(is.na(result$chromosome_count[2]))
  expect_equal(result$unfixable_error[2], 1L)
})

test_that("karyo_preprocessed: on_issues='stop' warning message cites count", {
  pp <- suppressMessages(preprocess_karyo(c("46,XX", NA, NA)))
  expect_warning(
    parse_karyo(pp),
    regexp = "2 of 3"
  )
})

test_that("karyo_preprocessed: clean-only input with on_issues='stop' is silent", {
  pp <- suppressMessages(preprocess_karyo(c("46,XX", "47,XY,+21")))
  expect_silent(parse_karyo(pp))
})

test_that("full pipeline: original_karyotype in parse output is the raw string, not preprocessed", {
  df <- data.frame(
    karyotype = c(".46,XX[20]"),
    stringsAsFactors = FALSE
  )
  ck <- check_karyo(df)
  pp <- suppressMessages(preprocess_karyo(ck))
  result <- suppressMessages(parse_karyo(pp))
  expect_equal(result$original_karyotype, ".46,XX[20]")
  expect_equal(result$preprocessed_karyotype, "46,XX[20]")
})

test_that("type guard: factor karyotype column triggers warning, still parses", {
  df <- data.frame(
    karyotype = factor(c("46,XX", "47,XY,+21")),
    stringsAsFactors = TRUE
  )
  expect_warning(
    result <- check_karyo(df),
    regexp = "not character"
  )
  expect_equal(nrow(result), 2L)
})

test_that("type guard: integer column named karyotype triggers warning", {
  df <- data.frame(karyotype = c(1L, 2L, 3L))
  expect_warning(
    check_karyo(df),
    regexp = "not character"
  )
})

test_that("type guard: factor karyotype column in parse_karyo triggers warning", {
  df <- data.frame(
    karyotype = factor(c("46,XX", "47,XY,+21")),
    stringsAsFactors = TRUE
  )
  expect_warning(
    result <- suppressMessages(parse_karyo(df)),
    regexp = "not character"
  )
  expect_equal(nrow(result), 2L)
})

test_that("parse_karyo: karyo_preprocessed accepted without karyotype_column", {
  pp <- suppressMessages(preprocess_karyo(c("46,XX", "47,XY,+21")))
  result <- suppressMessages(parse_karyo(pp))
  expect_equal(nrow(result), 2L)
  expect_false(is.na(result$chromosome_count[1]))
  expect_false(is.na(result$chromosome_count[2]))
})

test_that("parse_karyo: karyo_preprocessed propagates id without re-specification", {
  df <- data.frame(
    sample_id = c("A", "B"),
    karyotype = c("46,XX", "47,XY,+21"),
    stringsAsFactors = FALSE
  )
  pp <- suppressMessages(preprocess_karyo(check_karyo(df)))
  result <- suppressMessages(parse_karyo(pp))
  expect_true("sample_id" %in% names(result))
  expect_equal(result$sample_id, c("A", "B"))
})

test_that("parse_karyo: original_karyotype shows pre-fix string when using karyo_preprocessed", {
  raw <- ".46,XX,t(9;22)(q34;q11)[20]"
  pp <- suppressMessages(preprocess_karyo(raw))
  result <- suppressMessages(parse_karyo(pp))
  expect_equal(result$original_karyotype, raw)
  expect_equal(result$preprocessed_karyotype, "46,XX,t(9;22)(q34;q11)[20]")
})

test_that("check_karyo: verbose=FALSE produces no messages", {
  expect_silent(check_karyo(c("46,XX", "47,XY,+21"), verbose = FALSE))
})

test_that("check_karyo: verbose=TRUE produces messages and returns visibly", {
  msgs <- character(0)
  withCallingHandlers(
    {
      result <- check_karyo(c("46,XX", ".46,XY"), verbose = TRUE)
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_true(length(msgs) > 0)
  expect_s3_class(result, "karyo_check")
})

test_that("preprocess_karyo: verbose=FALSE produces no messages", {
  expect_silent(preprocess_karyo(c("46,XX", "47,XY,+21"), verbose = FALSE))
})

test_that("preprocess_karyo: verbose=TRUE produces messages", {
  suppressMessages(expect_message(
    preprocess_karyo(c("46,XX", ".46,XY"), verbose = TRUE),
    regexp = "Preprocessing"
  ))
})
