test_that("t(15;17)(q24;q21) detected", {
  r <- pk("46,XX,t(15;17)(q24;q21)")
  expect_equal(r$t_15_17_q24_q21, 1L)
})

test_that("t(15;17) reversed chromosome order", {
  r <- pk("46,XX,t(17;15)(q21;q24)")
  expect_equal(r$t_15_17_q24_q21, 1L)
})

test_that("t(15;17) with q22 variant band", {
  r <- pk("46,XX,t(15;17)(q22;q21)")
  expect_equal(r$t_15_17_q24_q21, 1L)
})

test_that("t(15;17) negative case", {
  r <- pk("46,XX")
  expect_equal(r$t_15_17_q24_q21, 0L)
})

test_that("t(8;21)(q22;q22) detected and reversed", {
  r1 <- pk("46,XY,t(8;21)(q22;q22)")
  expect_equal(r1$t_8_21_q22_q22, 1L)
  r2 <- pk("46,XY,t(21;8)(q22;q22)")
  expect_equal(r2$t_8_21_q22_q22, 1L)
})

test_that("t(8;21) with q21.3 sub-band variant", {
  r <- pk("46,XY,t(8;21)(q21.3;q22)")
  expect_equal(r$t_8_21_q22_q22, 1L)
})

test_that("inv(16)(p13q22) detected", {
  r <- pk("46,XX,inv(16)(p13q22)")
  expect_equal(r$inv_16_p13q22, 1L)
})

test_that("t(16;16)(p13;q22) detected", {
  r <- pk("46,XX,t(16;16)(p13;q22)")
  expect_equal(r$t_16_16_p13_q22, 1L)
})

test_that("t(9;11)(p21;q23) detected and reversed", {
  r1 <- pk("46,XX,t(9;11)(p21;q23)")
  expect_equal(r1$t_9_11_p21_q23, 1L)
  r2 <- pk("46,XX,t(11;9)(q23;p21)")
  expect_equal(r2$t_9_11_p21_q23, 1L)
})

test_that("t(6;9)(p22;q34) detected and reversed", {
  r1 <- pk("46,XY,t(6;9)(p22;q34)")
  expect_equal(r1$t_6_9_p22_q34, 1L)
  r2 <- pk("46,XY,t(9;6)(q34;p22)")
  expect_equal(r2$t_6_9_p22_q34, 1L)
})

test_that("inv(3)(q21q26) detected", {
  r <- pk("46,XY,inv(3)(q21q26)")
  expect_equal(r$inv_3_q21q26, 1L)
})

test_that("t(3;3)(q21;q26) detected", {
  r <- pk("46,XY,t(3;3)(q21;q26)")
  expect_equal(r$t_3_3_q21_q26, 1L)
})

test_that("t(9;22)(q34;q11) detected and reversed", {
  r1 <- pk("46,XY,t(9;22)(q34;q11)")
  expect_equal(r1$t_9_22_q34_q11, 1L)
  r2 <- pk("46,XY,t(22;9)(q11;q34)")
  expect_equal(r2$t_9_22_q34_q11, 1L)
})

test_that("batch: remaining specific translocations", {
  cases <- list(
    list(k = "46,XY,t(1;3)(p36;q21)", flag = "t_1_3_p36_q21"),
    list(k = "46,XX,t(1;22)(p13;q13)", flag = "t_1_22_p13_q13"),
    list(k = "46,XY,t(3;5)(q25;q35)", flag = "t_3_5_q25_q35"),
    list(k = "46,XX,t(5;11)(q35;p15)", flag = "t_5_11_q35_p15"),
    list(k = "46,XY,t(7;12)(q36;p13)", flag = "t_7_12_q36_p13"),
    list(k = "46,XX,t(8;16)(p11;p13)", flag = "t_8_16_p11_p13"),
    list(k = "46,XY,t(10;11)(p12;q14)", flag = "t_10_11_p12_q14"),
    list(k = "46,XX,t(11;12)(p15;p13)", flag = "t_11_12_p15_p13"),
    list(k = "46,XY,t(16;21)(p11;q22)", flag = "t_16_21_p11_q22"),
    list(k = "46,XX,t(16;21)(q24;q22)", flag = "t_16_21_q24_q22"),
    list(k = "46,XY,inv(16)(p13q24)", flag = "inv_16_p13q24")
  )
  for (tc in cases) {
    r <- pk(tc$k)
    expect_equal(r[[tc$flag]], 1L, info = paste("Testing", tc$flag))
  }
})

test_that("t(6;9)_other with non-canonical bands", {
  r <- pk("46,XY,t(6;9)(p23;q33)")
  expect_equal(r$t_6_9_other, 1L)
  expect_equal(r$t_6_9_p22_q34, 0L)
})

test_that("t(9;11)_other with non-canonical bands", {
  r <- pk("46,XX,t(9;11)(p22;q25)")
  expect_equal(r$t_9_11_other, 1L)
  expect_equal(r$t_9_11_p21_q23, 0L)
})

test_that("t(9;22)_other with non-canonical bands", {
  r <- pk("46,XY,t(9;22)(q32;q12)")
  expect_equal(r$t_9_22_other, 1L)
  expect_equal(r$t_9_22_q34_q11, 0L)
})

test_that("inv(3)_other with non-canonical bands", {
  r <- pk("46,XY,inv(3)(q22q28)")
  expect_equal(r$inv_3_other, 1L)
  expect_equal(r$inv_3_q21q26, 0L)
})

test_that("t(3;3)_other with non-canonical bands", {
  r <- pk("46,XY,t(3;3)(q25;q29)")
  expect_equal(r$t_3_3_other, 1L)
  expect_equal(r$t_3_3_q21_q26, 0L)
})

test_that("t(v;11p15) with arbitrary partner", {
  r <- pk("46,XX,t(4;11)(q21;p15)")
  expect_equal(r$t_v_11p15, 1L)
})

test_that("t(v;11q23) with arbitrary partner", {
  r <- pk("46,XY,t(6;11)(q27;q23)")
  expect_equal(r$t_v_11q23, 1L)
})

test_that("t(3q26;v) with arbitrary partner", {
  r <- pk("46,XX,t(3;7)(q26;p22)")
  expect_equal(r$t_3q26_v, 1L)
})

test_that("t(v;11q23) reversed orientation", {
  r <- pk("46,XX,t(11;4)(q23;q21)")
  expect_equal(r$t_v_11q23, 1L)
})

test_that("del(5q) detected", {
  r <- pk("46,XX,del(5)(q13q33)")
  expect_equal(r$del_5q, 1L)
})

test_that("del(5q) shorthand notation detected", {
  r <- pk("46,XY,del(5q)")
  expect_equal(r$del_5q, 1L)
})

test_that("del(5)(p14) does NOT trigger del(5q)", {
  r <- pk("46,XX,del(5)(p14)")
  expect_equal(r$del_5q, 0L)
})

test_that("add(5q) detected", {
  r <- pk("46,XY,add(5)(q31)")
  expect_equal(r$add_5q, 1L)
})

test_that("add(5q) shorthand notation detected", {
  r <- pk("46,XY,add(5q)")
  expect_equal(r$add_5q, 1L)
})

test_that("t(5q) detected", {
  r <- pk("46,XX,t(5;17)(q33;p13)")
  expect_equal(r$t_5q, 1L)
})

test_that("del(7q) detected", {
  r <- pk("46,XY,del(7)(q22)")
  expect_equal(r$del_7q, 1L)
})

test_that("del(7q) shorthand notation detected", {
  r <- pk("46,XY,del(7q)")
  expect_equal(r$del_7q, 1L)
})

test_that("del(12p) detected", {
  r <- pk("46,XX,del(12)(p12)")
  expect_equal(r$del_12p, 1L)
})

test_that("del(12p) shorthand notation detected", {
  r <- pk("46,XX,del(12p)")
  expect_equal(r$del_12p, 1L)
})

test_that("t(12p) detected", {
  r <- pk("46,XY,t(12;14)(p13;q32)")
  expect_equal(r$t_12p, 1L)
})

test_that("add(12p) detected", {
  r <- pk("46,XX,add(12)(p11)")
  expect_equal(r$add_12p, 1L)
})

test_that("add(12p) shorthand notation detected", {
  r <- pk("46,XX,add(12p)")
  expect_equal(r$add_12p, 1L)
})

test_that("del(13q) detected", {
  r <- pk("46,XY,del(13)(q14)")
  expect_equal(r$del_13q, 1L)
})

test_that("del(13q) shorthand notation detected", {
  r <- pk("46,XY,del(13q)")
  expect_equal(r$del_13q, 1L)
})

test_that("i(17q) detected", {
  r <- pk("46,XX,i(17)(q10)")
  expect_equal(r$i_17q, 1L)
})

test_that("i(17q) shorthand notation detected", {
  r <- pk("46,XX,i(17q)")
  expect_equal(r$i_17q, 1L)
})

test_that("add(17p) detected", {
  r <- pk("46,XY,add(17)(p11)")
  expect_equal(r$add_17p, 1L)
})

test_that("add(17p) shorthand notation detected", {
  r <- pk("46,XY,add(17p)")
  expect_equal(r$add_17p, 1L)
})

test_that("del(17p) detected", {
  r <- pk("46,XX,del(17)(p13)")
  expect_equal(r$del_17p, 1L)
})

test_that("del(17p) shorthand notation detected", {
  r <- pk("46,XX,del(17p)")
  expect_equal(r$del_17p, 1L)
})

test_that("del(20q) detected", {
  r <- pk("46,XY,del(20)(q11)")
  expect_equal(r$del_20q, 1L)
})

test_that("del(20q) shorthand notation detected", {
  r <- pk("46,XY,del(20q)")
  expect_equal(r$del_20q, 1L)
})

test_that("del(11q) detected", {
  r <- pk("46,XX,del(11)(q23)")
  expect_equal(r$del_11q, 1L)
})

test_that("del(11q) shorthand notation detected", {
  r <- pk("46,XX,del(11q)")
  expect_equal(r$del_11q, 1L)
})

test_that("idic(X)(q13) detected", {
  r <- pk("46,XX,idic(X)(q13)")
  expect_equal(r$idic_X_q13, 1L)
})

test_that("idic(X)(q13.1) sub-band still matches idic(X)(q13) via prefix", {
  r <- pk("46,XX,idic(X)(q13.1)")
  expect_equal(r$idic_X_q13, 1L)
  expect_equal(r$general_isodicentric, 1L)
})

test_that("idic(X)(q14) does NOT match idic(X)(q13); fires isodicentric instead", {
  r <- pk("46,XX,idic(X)(q14)")
  expect_equal(r$idic_X_q13, 0L)
  expect_equal(r$general_isodicentric, 1L)
})

test_that("del(7)(p11) does NOT trigger del(7q)", {
  r <- pk("46,XY,del(7)(p11)")
  expect_equal(r$del_7q, 0L)
})

test_that("del(12)(q22) does NOT trigger del(12p)", {
  r <- pk("46,XX,del(12)(q22)")
  expect_equal(r$del_12p, 0L)
})

test_that("dicentric detected", {
  r <- pk("46,XY,dic(7;9)(p11;p13)")
  expect_equal(r$general_dicentric, 1L)
})

test_that("isodicentric detected", {
  r <- pk("46,XX,idic(7)(p11)")
  expect_equal(r$general_isodicentric, 1L)
})

test_that("pseudodicentric detected", {
  r <- pk("46,XX,psu dic(15;22)(q11;p11)")
  expect_equal(r$general_pseudodicentric, 1L)
})

test_that("ring chromosome detected", {
  r <- pk("46,XX,r(7)(p22q36)")
  expect_equal(r$general_ring, 1L)
})

test_that("insertion detected", {
  r <- pk("46,XY,ins(5;11)(p14;q13q23)")
  expect_equal(r$general_insertion, 1L)
})

test_that("duplication detected", {
  r <- pk("46,XX,dup(1)(q21q32)")
  expect_equal(r$general_duplication, 1L)
})

test_that("triplication detected", {
  r <- pk("46,XX,trp(1)(q21q32)")
  expect_equal(r$general_triplication, 1L)
})

test_that("general_addition detected", {
  r <- pk("46,XX,add(1)(p36)")
  expect_equal(r$general_addition, 1L)
})

test_that("general_inversion detected", {
  r <- pk("46,XY,inv(9)(p11q13)")
  expect_equal(r$general_inversion, 1L)
})

test_that("general_deletion detected", {
  r <- pk("46,XX,del(1)(q21)")
  expect_equal(r$general_deletion, 1L)
})

test_that("marker_chromosome detected", {
  r <- pk("47,XY,+mar")
  expect_equal(r$general_marker, 1L)
})

test_that("general_translocation detected", {
  r <- pk("46,XX,t(2;7)(p11;q22)")
  expect_equal(r$general_translocation, 1L)
})

test_that("general_translocation negative", {
  r <- pk("46,XX")
  expect_equal(r$general_translocation, 0L)
})

test_that("general_translocation: three-way translocation detected", {
  r <- pk("46,XX,t(1;2;3)(p11;q22;p13)")
  expect_equal(r$general_translocation, 1L)
})

test_that("general_translocation: X chromosome partner detected", {
  r <- pk("46,X,t(X;1)(p22;q11)")
  expect_equal(r$general_translocation, 1L)
})

test_that("general_translocation: no band info detected", {
  r <- pk("46,XX,t(1;2)")
  expect_equal(r$general_translocation, 1L)
})

test_that("general_translocation: der with embedded t does NOT match", {
  r <- pk("46,XX,der(1)t(1;2)(p11;q22)")
  expect_equal(r$general_translocation, 0L)
})

test_that("general_translocation: Y chromosome as first partner detected", {
  r <- pk("46,XY,t(Y;1)(q11;p11)")
  expect_equal(r$general_translocation, 1L)
})

test_that("specific inv(3) and general_inversion co-fire", {
  r <- pk("46,XY,inv(3)(q21q26)")
  expect_equal(r$inv_3_q21q26, 1L)
  expect_equal(r$general_inversion, 1L)
})

test_that("specific del(5q) and general_deletion co-fire", {
  r <- pk("46,XX,del(5)(q13q33)")
  expect_equal(r$del_5q, 1L)
  expect_equal(r$general_deletion, 1L)
})

test_that("idic(X)(q13) and general_isodicentric co-fire", {
  r <- pk("46,XX,idic(X)(q13)")
  expect_equal(r$idic_X_q13, 1L)
  expect_equal(r$general_isodicentric, 1L)
})

test_that("canonical t(9;22) fires, _other variant does NOT", {
  r <- pk("46,XY,t(9;22)(q34;q11)")
  expect_equal(r$t_9_22_q34_q11, 1L)
  expect_equal(r$t_9_22_other, 0L)
})

test_that("_other excludes its canonical breakpoints but fires for the rest", {
  # mutual exclusivity is in the regex (negative lookahead), not priority
  nc <- pk("46,XY,t(9;22)(q11;q34)")
  expect_equal(nc$t_9_22_q34_q11, 0L)
  expect_equal(nc$t_9_22_other, 1L)
  expect_equal(pk("46,XY,inv(3)(q24q26)")$inv_3_other, 1L)
  expect_equal(pk("46,XY,inv(3)(q21q26)")$inv_3_other, 0L)
})

test_that("family flags co-fire with the specific rule on the same token", {
  r1 <- pk("46,XX,t(9;11)(p21;q23)")
  expect_equal(r1$t_9_11_p21_q23, 1L)
  expect_equal(r1$t_v_11q23, 1L)
  r2 <- pk("46,XX,t(3;5)(q25;q35)")
  expect_equal(r2$t_3_5_q25_q35, 1L)
  expect_equal(r2$t_5q, 1L)
})

test_that("specific tx and general_translocation co-fire", {
  r <- pk("46,XY,t(9;22)(q34;q11)")
  expect_equal(r$t_9_22_q34_q11, 1L)
  expect_equal(r$general_translocation, 1L)
})

test_that("variable partner and general_translocation co-fire", {
  r <- pk("46,XY,t(6;11)(q27;q23)")
  expect_equal(r$t_v_11q23, 1L)
  expect_equal(r$general_translocation, 1L)
})

test_that("chromosome-specific tx and general_translocation co-fire", {
  r <- pk("46,XX,t(5;17)(q33;p13)")
  expect_equal(r$t_5q, 1L)
  expect_equal(r$general_translocation, 1L)
})

test_that("t(v;11p15) and general_translocation co-fire", {
  r <- pk("46,XX,t(4;11)(q21;p15)")
  expect_equal(r$t_v_11p15, 1L)
  expect_equal(r$general_translocation, 1L)
})

test_that("der with embedded del(5q): both del(5q) and general_derivative fire", {
  r <- pk("46,XX,der(5)del(5)(q11q34)")
  expect_equal(r$del_5q, 1L)
  expect_equal(r$general_derivative, 1L)
})

test_that("der of other chrom with embedded del(5q): both fire", {
  r <- pk("46,XX,der(3)del(5)(q11q34)")
  expect_equal(r$del_5q, 1L)
  expect_equal(r$general_derivative, 1L)
})

test_that("der with embedded ins and del: insertion, del(5q), and general_derivative all fire", {
  r <- pk("46,XX,der(5)ins(5;17)(p11;??)del(5)(q11)")
  expect_equal(r$general_insertion, 1L)
  expect_equal(r$del_5q, 1L)
  expect_equal(r$general_derivative, 1L)
})

test_that("der with unspecific embedded del: general_deletion and general_derivative fire", {
  r <- pk("46,XX,der(9)del(9)(p11)")
  expect_equal(r$general_deletion, 1L)
  expect_equal(r$general_derivative, 1L)
})

test_that("der(8)t(8;21) with monosomy: t(8;21) and general_derivative fire, monosomal overridden to 0", {
  r <- pk("45,XX,-7,der(8)t(8;21)(q22;q22)")
  expect_equal(r$t_8_21_q22_q22, 1L)
  expect_equal(r$general_derivative, 1L)
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("der(8)t(8;21): general_translocation stays 0 (der token, not bare t)", {
  r <- pk("45,XX,-7,der(8)t(8;21)(q22;q22)")
  expect_equal(r$general_translocation, 0L)
})

test_that("standalone t(9;22): specific and general_translocation co-fire", {
  r <- pk("46,XX,t(9;22)(q34;q11)")
  expect_equal(r$t_9_22_q34_q11, 1L)
  expect_equal(r$general_translocation, 1L)
})

test_that("myeloid_rules has karyo_rules class and expected columns", {
  expect_s3_class(myeloid_rules, "karyo_rules")
  expect_true(tibble::is_tibble(myeloid_rules))
  expect_true(nrow(myeloid_rules) > 0L)
  expect_true(all(
    c(
      "flag_name",
      "regex"
    ) %in%
      names(myeloid_rules)
  ))
})

test_that("a punctuation-heavy flag_name maps to a sanitized, backtick-free column name", {
  custom <- validate_rules(data.frame(
    flag_name = "t(9;22)(q34;q11)",
    regex = "t\\(9;22\\)\\(q34;q11\\)"
  ))
  r <- parse_karyo(
    "46,XX,t(9;22)(q34;q11)",
    rules = custom,
    on_issues = "warn",
    verbose = FALSE
  )
  expect_true("t_9_22_q34_q11" %in% names(r))
  expect_equal(r$t_9_22_q34_q11, 1L)
})

test_that("general structural flags fire independently of the rule set", {
  # general_* are detected by the parser, not by myeloid_rules
  expect_equal(pk("46,XX,dic(1;7)(p11;p11)")$general_dicentric, 1L)
  expect_equal(pk("47,XY,+mar")$general_marker, 1L)
})

test_that("general_isochromosome fires for any isochromosome and co-fires with i(17q)", {
  expect_equal(pk("46,XX,i(7)(q10)")$general_isochromosome, 1L)
  r <- pk("46,XX,i(17)(q10)")
  expect_equal(r$i_17q, 1L)
  expect_equal(r$general_isochromosome, 1L)
  # isodicentric tokens are not isochromosomes
  expect_equal(pk("46,XX,idic(X)(q13)")$general_isochromosome, 0L)
})

test_that("general_* still fire with a custom rule set lacking general rules", {
  one_rule <- validate_rules(data.frame(
    flag_name = "t(9;22)(q34;q11)",
    regex = "t\\(9;22\\)\\(q34;q11\\)",
    priority = 100,
    competition_group = "translocation"
  ))
  r <- parse_karyo(
    "46,XX,t(9;22)(q34;q11),del(1)(q21)",
    rules = one_rule,
    on_issues = "warn",
    verbose = FALSE
  )
  expect_equal(r$t_9_22_q34_q11, 1L)
  expect_equal(r$general_translocation, 1L)
  expect_equal(r$general_deletion, 1L)
})

test_that("validate_rules() attaches karyo_rules class to a valid data frame", {
  df <- as.data.frame(myeloid_rules)
  class(df) <- setdiff(class(df), "karyo_rules")
  vr <- validate_rules(df)
  expect_s3_class(vr, "karyo_rules")
})

test_that("validate_rules() errors on missing columns", {
  bad <- myeloid_rules[, "flag_name"]
  class(bad) <- setdiff(class(bad), "karyo_rules")
  expect_error(validate_rules(bad), "missing required columns")
})

test_that("validate_rules() errors on non-data-frame input", {
  expect_error(validate_rules("not a data frame"), "must be a data frame")
})

test_that("parse_karyo() errors when rules is a plain data frame", {
  df <- as.data.frame(myeloid_rules)
  class(df) <- setdiff(class(df), "karyo_rules")
  expect_error(parse_karyo("46,XX", rules = df), "karyo_rules")
})

test_that("rules accepts a list of karyo_rules tables and combines their flags", {
  set_a <- validate_rules(data.frame(
    flag_name = "flagA",
    regex = "t\\(9;22\\)\\(q34;q11\\)"
  ))
  set_b <- validate_rules(data.frame(
    flag_name = "flagB",
    regex = "del\\(5\\)\\(q13q33\\)"
  ))
  r <- parse_karyo(
    c("46,XX,t(9;22)(q34;q11)", "46,XX,del(5)(q13q33)"),
    rules = list(set_a, set_b),
    on_issues = "warn",
    verbose = FALSE
  )
  expect_true(all(c("flagA", "flagB") %in% names(r)))
  expect_equal(r$flagA, c(1L, 0L))
  expect_equal(r$flagB, c(0L, 1L))
})

test_that("a single-element rules list behaves like passing the table directly", {
  r_list <- parse_karyo(
    "46,XX,t(9;22)(q34;q11)",
    rules = list(myeloid_rules),
    on_issues = "warn",
    verbose = FALSE
  )
  r_direct <- pk("46,XX,t(9;22)(q34;q11)")
  expect_identical(r_list, r_direct)
})

test_that("rules list: a flag_name shared across tables errors instead of merging", {
  set_a <- validate_rules(data.frame(
    flag_name = "shared_flag",
    regex = "t\\(9;22\\)\\(q34;q11\\)"
  ))
  set_b <- validate_rules(data.frame(
    flag_name = "shared_flag",
    regex = "del\\(5\\)\\(q13q33\\)"
  ))
  expect_error(
    parse_karyo(
      c("46,XX,t(9;22)(q34;q11)", "46,XX,del(5)(q13q33)", "46,XX"),
      rules = list(set_a, set_b),
      on_issues = "warn",
      verbose = FALSE
    ),
    "duplicate `flag_name`"
  )
})

test_that("validate_rules() errors on duplicate flag_name within a single table", {
  bad <- data.frame(
    flag_name = c("del(5q)", "del(5q)"),
    regex = c("del\\(5q", "del\\(5\\)\\(q")
  )
  expect_error(validate_rules(bad), "duplicate `flag_name`")
})

test_that("validate_rules() errors when distinct flag_names collide once sanitized", {
  bad <- data.frame(
    flag_name = c("t(9;22)", "t(9,22)"),
    regex = c("t\\(9;22\\)", "t\\(9;22\\)")
  )
  expect_error(validate_rules(bad), "collide once sanitized")
})

test_that("rules errors on an empty list", {
  expect_error(
    parse_karyo("46,XX", rules = list()),
    "at least one"
  )
})

test_that("rules errors when a list element is not a karyo_rules object", {
  expect_error(
    parse_karyo("46,XX", rules = list(myeloid_rules, "not a karyo_rules")),
    "Every element of `rules`"
  )
})

test_that("preprocessed_karyotype under on_issues='warn' is normalize_iscn() of input", {
  r <- parse_karyo(c("46,XX", "46,XY,+8"), on_issues = "warn", verbose = FALSE)
  expect_equal(r$preprocessed_karyotype[1], "46,XX")
  expect_equal(r$preprocessed_karyotype[2], "46,XY,+8")
})

test_that("preprocessed_karyotype under on_issues='fix' equals the cleaned string", {
  r <- parse_karyo(".46,XX", on_issues = "preprocess", verbose = FALSE)
  expect_equal(r$preprocessed_karyotype, "46,XX")
})

test_that("preprocessed_karyotype column is present in output", {
  r <- pk("46,XX")
  expect_true("preprocessed_karyotype" %in% names(r))
})
