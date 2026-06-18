test_that("t(15;17)(q24;q21) detected", {
  r <- pk("46,XX,t(15;17)(q24;q21)")
  expect_equal(r$`t(15;17)(q24;q21)`, 1L)
})

test_that("t(15;17) reversed chromosome order", {
  r <- pk("46,XX,t(17;15)(q21;q24)")
  expect_equal(r$`t(15;17)(q24;q21)`, 1L)
})

test_that("t(15;17) with q22 variant band", {
  r <- pk("46,XX,t(15;17)(q22;q21)")
  expect_equal(r$`t(15;17)(q24;q21)`, 1L)
})

test_that("t(15;17) negative case", {
  r <- pk("46,XX")
  expect_equal(r$`t(15;17)(q24;q21)`, 0L)
})

test_that("t(8;21)(q22;q22) detected and reversed", {
  r1 <- pk("46,XY,t(8;21)(q22;q22)")
  expect_equal(r1$`t(8;21)(q22;q22)`, 1L)
  r2 <- pk("46,XY,t(21;8)(q22;q22)")
  expect_equal(r2$`t(8;21)(q22;q22)`, 1L)
})

test_that("t(8;21) with q21.3 sub-band variant", {
  r <- pk("46,XY,t(8;21)(q21.3;q22)")
  expect_equal(r$`t(8;21)(q22;q22)`, 1L)
})

test_that("inv(16)(p13q22) detected", {
  r <- pk("46,XX,inv(16)(p13q22)")
  expect_equal(r$`inv(16)(p13q22)`, 1L)
})

test_that("t(16;16)(p13;q22) detected", {
  r <- pk("46,XX,t(16;16)(p13;q22)")
  expect_equal(r$`t(16;16)(p13;q22)`, 1L)
})

test_that("t(9;11)(p21;q23) detected and reversed", {
  r1 <- pk("46,XX,t(9;11)(p21;q23)")
  expect_equal(r1$`t(9;11)(p21;q23)`, 1L)
  r2 <- pk("46,XX,t(11;9)(q23;p21)")
  expect_equal(r2$`t(9;11)(p21;q23)`, 1L)
})

test_that("t(6;9)(p22;q34) detected and reversed", {
  r1 <- pk("46,XY,t(6;9)(p22;q34)")
  expect_equal(r1$`t(6;9)(p22;q34)`, 1L)
  r2 <- pk("46,XY,t(9;6)(q34;p22)")
  expect_equal(r2$`t(6;9)(p22;q34)`, 1L)
})

test_that("inv(3)(q21q26) detected", {
  r <- pk("46,XY,inv(3)(q21q26)")
  expect_equal(r$`inv(3)(q21q26)`, 1L)
})

test_that("t(3;3)(q21;q26) detected", {
  r <- pk("46,XY,t(3;3)(q21;q26)")
  expect_equal(r$`t(3;3)(q21;q26)`, 1L)
})

test_that("t(9;22)(q34;q11) detected and reversed", {
  r1 <- pk("46,XY,t(9;22)(q34;q11)")
  expect_equal(r1$`t(9;22)(q34;q11)`, 1L)
  r2 <- pk("46,XY,t(22;9)(q11;q34)")
  expect_equal(r2$`t(9;22)(q34;q11)`, 1L)
})

test_that("batch: remaining specific translocations", {
  cases <- list(
    list(k = "46,XY,t(1;3)(p36;q21)", flag = "t(1;3)(p36;q21)"),
    list(k = "46,XX,t(1;22)(p13;q13)", flag = "t(1;22)(p13;q13)"),
    list(k = "46,XY,t(3;5)(q25;q35)", flag = "t(3;5)(q25;q35)"),
    list(k = "46,XX,t(5;11)(q35;p15)", flag = "t(5;11)(q35;p15)"),
    list(k = "46,XY,t(7;12)(q36;p13)", flag = "t(7;12)(q36;p13)"),
    list(k = "46,XX,t(8;16)(p11;p13)", flag = "t(8;16)(p11;p13)"),
    list(k = "46,XY,t(10;11)(p12;q14)", flag = "t(10;11)(p12;q14)"),
    list(k = "46,XX,t(11;12)(p15;p13)", flag = "t(11;12)(p15;p13)"),
    list(k = "46,XY,t(16;21)(p11;q22)", flag = "t(16;21)(p11;q22)"),
    list(k = "46,XX,t(16;21)(q24;q22)", flag = "t(16;21)(q24;q22)"),
    list(k = "46,XY,inv(16)(p13q24)", flag = "inv(16)(p13q24)")
  )
  for (tc in cases) {
    r <- pk(tc$k)
    expect_equal(r[[tc$flag]], 1L, info = paste("Testing", tc$flag))
  }
})

test_that("t(6;9)_other with non-canonical bands", {
  r <- pk("46,XY,t(6;9)(p23;q33)")
  expect_equal(r$`t(6;9)_other`, 1L)
  expect_equal(r$`t(6;9)(p22;q34)`, 0L)
})

test_that("t(9;11)_other with non-canonical bands", {
  r <- pk("46,XX,t(9;11)(p22;q25)")
  expect_equal(r$`t(9;11)_other`, 1L)
  expect_equal(r$`t(9;11)(p21;q23)`, 0L)
})

test_that("t(9;22)_other with non-canonical bands", {
  r <- pk("46,XY,t(9;22)(q32;q12)")
  expect_equal(r$`t(9;22)_other`, 1L)
  expect_equal(r$`t(9;22)(q34;q11)`, 0L)
})

test_that("inv(3)_other with non-canonical bands", {
  r <- pk("46,XY,inv(3)(q22q28)")
  expect_equal(r$`inv(3)_other`, 1L)
  expect_equal(r$`inv(3)(q21q26)`, 0L)
})

test_that("t(3;3)_other with non-canonical bands", {
  r <- pk("46,XY,t(3;3)(q25;q29)")
  expect_equal(r$`t(3;3)_other`, 1L)
  expect_equal(r$`t(3;3)(q21;q26)`, 0L)
})

test_that("t(v;11p15) with arbitrary partner", {
  r <- pk("46,XX,t(4;11)(q21;p15)")
  expect_equal(r$`t(v;11p15)`, 1L)
})

test_that("t(v;11q23) with arbitrary partner", {
  r <- pk("46,XY,t(6;11)(q27;q23)")
  expect_equal(r$`t(v;11q23)`, 1L)
})

test_that("t(3q26;v) with arbitrary partner", {
  r <- pk("46,XX,t(3;7)(q26;p22)")
  expect_equal(r$`t(3q26;v)`, 1L)
})

test_that("t(v;11q23) reversed orientation", {
  r <- pk("46,XX,t(11;4)(q23;q21)")
  expect_equal(r$`t(v;11q23)`, 1L)
})

test_that("del(5q) detected", {
  r <- pk("46,XX,del(5)(q13q33)")
  expect_equal(r$`del(5q)`, 1L)
})

test_that("del(5q) shorthand notation detected", {
  r <- pk("46,XY,del(5q)")
  expect_equal(r$`del(5q)`, 1L)
})

test_that("del(5)(p14) does NOT trigger del(5q)", {
  r <- pk("46,XX,del(5)(p14)")
  expect_equal(r$`del(5q)`, 0L)
})

test_that("add(5q) detected", {
  r <- pk("46,XY,add(5)(q31)")
  expect_equal(r$`add(5q)`, 1L)
})

test_that("add(5q) shorthand notation detected", {
  r <- pk("46,XY,add(5q)")
  expect_equal(r$`add(5q)`, 1L)
})

test_that("t(5q) detected", {
  r <- pk("46,XX,t(5;17)(q33;p13)")
  expect_equal(r$`t(5q)`, 1L)
})

test_that("del(7q) detected", {
  r <- pk("46,XY,del(7)(q22)")
  expect_equal(r$`del(7q)`, 1L)
})

test_that("del(7q) shorthand notation detected", {
  r <- pk("46,XY,del(7q)")
  expect_equal(r$`del(7q)`, 1L)
})

test_that("del(12p) detected", {
  r <- pk("46,XX,del(12)(p12)")
  expect_equal(r$`del(12p)`, 1L)
})

test_that("del(12p) shorthand notation detected", {
  r <- pk("46,XX,del(12p)")
  expect_equal(r$`del(12p)`, 1L)
})

test_that("t(12p) detected", {
  r <- pk("46,XY,t(12;14)(p13;q32)")
  expect_equal(r$`t(12p)`, 1L)
})

test_that("add(12p) detected", {
  r <- pk("46,XX,add(12)(p11)")
  expect_equal(r$`add(12p)`, 1L)
})

test_that("add(12p) shorthand notation detected", {
  r <- pk("46,XX,add(12p)")
  expect_equal(r$`add(12p)`, 1L)
})

test_that("del(13q) detected", {
  r <- pk("46,XY,del(13)(q14)")
  expect_equal(r$`del(13q)`, 1L)
})

test_that("del(13q) shorthand notation detected", {
  r <- pk("46,XY,del(13q)")
  expect_equal(r$`del(13q)`, 1L)
})

test_that("i(17q) detected", {
  r <- pk("46,XX,i(17)(q10)")
  expect_equal(r$`i(17q)`, 1L)
})

test_that("i(17q) shorthand notation detected", {
  r <- pk("46,XX,i(17q)")
  expect_equal(r$`i(17q)`, 1L)
})

test_that("add(17p) detected", {
  r <- pk("46,XY,add(17)(p11)")
  expect_equal(r$`add(17p)`, 1L)
})

test_that("add(17p) shorthand notation detected", {
  r <- pk("46,XY,add(17p)")
  expect_equal(r$`add(17p)`, 1L)
})

test_that("del(17p) detected", {
  r <- pk("46,XX,del(17)(p13)")
  expect_equal(r$`del(17p)`, 1L)
})

test_that("del(17p) shorthand notation detected", {
  r <- pk("46,XX,del(17p)")
  expect_equal(r$`del(17p)`, 1L)
})

test_that("del(20q) detected", {
  r <- pk("46,XY,del(20)(q11)")
  expect_equal(r$`del(20q)`, 1L)
})

test_that("del(20q) shorthand notation detected", {
  r <- pk("46,XY,del(20q)")
  expect_equal(r$`del(20q)`, 1L)
})

test_that("del(11q) detected", {
  r <- pk("46,XX,del(11)(q23)")
  expect_equal(r$`del(11q)`, 1L)
})

test_that("del(11q) shorthand notation detected", {
  r <- pk("46,XX,del(11q)")
  expect_equal(r$`del(11q)`, 1L)
})

test_that("idic(X)(q13) detected", {
  r <- pk("46,XX,idic(X)(q13)")
  expect_equal(r$`idic(X)(q13)`, 1L)
})

test_that("idic(X)(q13.1) sub-band still matches idic(X)(q13) via prefix", {
  r <- pk("46,XX,idic(X)(q13.1)")
  expect_equal(r$`idic(X)(q13)`, 1L)
  expect_equal(r$isodicentric, 0L)
})

test_that("idic(X)(q14) does NOT match idic(X)(q13); fires isodicentric instead", {
  r <- pk("46,XX,idic(X)(q14)")
  expect_equal(r$`idic(X)(q13)`, 0L)
  expect_equal(r$isodicentric, 1L)
})

test_that("del(7)(p11) does NOT trigger del(7q)", {
  r <- pk("46,XY,del(7)(p11)")
  expect_equal(r$`del(7q)`, 0L)
})

test_that("del(12)(q22) does NOT trigger del(12p)", {
  r <- pk("46,XX,del(12)(q22)")
  expect_equal(r$`del(12p)`, 0L)
})

test_that("dicentric detected", {
  r <- pk("46,XY,dic(7;9)(p11;p13)")
  expect_equal(r$dicentric, 1L)
})

test_that("isodicentric detected", {
  r <- pk("46,XX,idic(7)(p11)")
  expect_equal(r$isodicentric, 1L)
})

test_that("pseudodicentric detected", {
  r <- pk("46,XX,psu dic(15;22)(q11;p11)")
  expect_equal(r$pseudodicentric, 1L)
})

test_that("ring chromosome detected", {
  r <- pk("46,XX,r(7)(p22q36)")
  expect_equal(r$ring_chromosome, 1L)
})

test_that("insertion detected", {
  r <- pk("46,XY,ins(5;11)(p14;q13q23)")
  expect_equal(r$insertion, 1L)
})

test_that("duplication detected", {
  r <- pk("46,XX,dup(1)(q21q32)")
  expect_equal(r$duplication, 1L)
})

test_that("triplication detected", {
  r <- pk("46,XX,trp(1)(q21q32)")
  expect_equal(r$triplication, 1L)
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
  expect_equal(r$marker_chromosome, 1L)
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

test_that("specific inv(3) fires, general_inversion does NOT", {
  r <- pk("46,XY,inv(3)(q21q26)")
  expect_equal(r$`inv(3)(q21q26)`, 1L)
  expect_equal(r$general_inversion, 0L)
})

test_that("specific del(5q) fires, general_deletion does NOT", {
  r <- pk("46,XX,del(5)(q13q33)")
  expect_equal(r$`del(5q)`, 1L)
  expect_equal(r$general_deletion, 0L)
})

test_that("idic(X)(q13) fires, general isodicentric does NOT", {
  r <- pk("46,XX,idic(X)(q13)")
  expect_equal(r$`idic(X)(q13)`, 1L)
  expect_equal(r$isodicentric, 0L)
})

test_that("canonical t(9;22) fires, _other variant does NOT", {
  r <- pk("46,XY,t(9;22)(q34;q11)")
  expect_equal(r$`t(9;22)(q34;q11)`, 1L)
  expect_equal(r$`t(9;22)_other`, 0L)
})

test_that("specific tx fires, general_translocation does NOT", {
  r <- pk("46,XY,t(9;22)(q34;q11)")
  expect_equal(r$`t(9;22)(q34;q11)`, 1L)
  expect_equal(r$general_translocation, 0L)
})

test_that("variable partner fires, general_translocation does NOT", {
  r <- pk("46,XY,t(6;11)(q27;q23)")
  expect_equal(r$`t(v;11q23)`, 1L)
  expect_equal(r$general_translocation, 0L)
})

test_that("chromosome-specific tx fires, general_translocation does NOT", {
  r <- pk("46,XX,t(5;17)(q33;p13)")
  expect_equal(r$`t(5q)`, 1L)
  expect_equal(r$general_translocation, 0L)
})

test_that("t(v;11p15) fires, general_translocation does NOT", {
  r <- pk("46,XX,t(4;11)(q21;p15)")
  expect_equal(r$`t(v;11p15)`, 1L)
  expect_equal(r$general_translocation, 0L)
})

test_that("der with embedded del(5q): both del(5q) and derivative_chromosome fire", {
  r <- pk("46,XX,der(5)del(5)(q11q34)")
  expect_equal(r$`del(5q)`, 1L)
  expect_equal(r$derivative_chromosome, 1L)
})

test_that("der of other chrom with embedded del(5q): both fire", {
  r <- pk("46,XX,der(3)del(5)(q11q34)")
  expect_equal(r$`del(5q)`, 1L)
  expect_equal(r$derivative_chromosome, 1L)
})

test_that("der with embedded ins and del: insertion, del(5q), and derivative_chromosome all fire", {
  r <- pk("46,XX,der(5)ins(5;17)(p11;??)del(5)(q11)")
  expect_equal(r$insertion, 1L)
  expect_equal(r$`del(5q)`, 1L)
  expect_equal(r$derivative_chromosome, 1L)
})

test_that("der with unspecific embedded del: general_deletion and derivative_chromosome fire", {
  r <- pk("46,XX,der(9)del(9)(p11)")
  expect_equal(r$general_deletion, 1L)
  expect_equal(r$derivative_chromosome, 1L)
})

test_that("der(8)t(8;21) with monosomy: t(8;21) and derivative_chromosome fire, monosomal overridden to 0", {
  r <- pk("45,XX,-7,der(8)t(8;21)(q22;q22)")
  expect_equal(r$`t(8;21)(q22;q22)`, 1L)
  expect_equal(r$derivative_chromosome, 1L)
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("der(8)t(8;21): general_translocation suppressed within translocation group", {
  r <- pk("45,XX,-7,der(8)t(8;21)(q22;q22)")
  expect_equal(r$general_translocation, 0L)
})

test_that("standalone t(9;22): specificity within translocation group preserved", {
  r <- pk("46,XX,t(9;22)(q34;q11)")
  expect_equal(r$`t(9;22)(q34;q11)`, 1L)
  expect_equal(r$general_translocation, 0L)
})

test_that("myeloid_rules has karyo_rules class and expected columns", {
  expect_s3_class(myeloid_rules, "karyo_rules")
  expect_true(tibble::is_tibble(myeloid_rules))
  expect_true(nrow(myeloid_rules) > 0L)
  expect_true(all(
    c(
      "flag_name",
      "regex",
      "category",
      "priority",
      "counts_for_monosomal",
      "competition_group"
    ) %in%
      names(myeloid_rules)
  ))
})

test_that("base_rules is a karyo_rules of only general rules", {
  expect_s3_class(base_rules, "karyo_rules")
  expect_true(all(base_rules$category == "general"))
  expect_true(nrow(base_rules) > 0L)
})

test_that("myeloid_rules is base_rules plus myeloid-specific rules", {
  # every base rule is carried into myeloid_rules unchanged
  base_in_myeloid <- merge(
    as.data.frame(base_rules),
    as.data.frame(myeloid_rules)
  )
  expect_equal(nrow(base_in_myeloid), nrow(base_rules))
  # myeloid_rules also carries non-general (disease-specific) rules
  expect_true(any(myeloid_rules$category != "general"))
})

test_that("a base (general) rule still fires under myeloid_rules", {
  expect_equal(pk("46,XX,dic(1;7)(p11;p11)")$dicentric, 1L)
})

test_that("myeloid_rules priorities are positive numeric", {
  expect_true(
    is.integer(myeloid_rules$priority) || is.numeric(myeloid_rules$priority)
  )
  expect_true(all(myeloid_rules$priority >= 1L))
})

test_that("validate_rules() attaches karyo_rules class to a valid data frame", {
  df <- as.data.frame(myeloid_rules)
  class(df) <- setdiff(class(df), "karyo_rules")
  vr <- validate_rules(df)
  expect_s3_class(vr, "karyo_rules")
})

test_that("validate_rules() errors on missing columns", {
  bad <- myeloid_rules[, c("flag_name", "regex")]
  class(bad) <- setdiff(class(bad), "karyo_rules")
  expect_error(validate_rules(bad), "missing required columns")
})

test_that("validate_rules() errors on non-data-frame input", {
  expect_error(validate_rules("not a data frame"), "must be a data frame")
})

test_that("validate_rules() errors on invalid priority", {
  bad <- as.data.frame(myeloid_rules)
  class(bad) <- setdiff(class(bad), "karyo_rules")
  bad$priority[1] <- -1
  expect_error(validate_rules(bad), "priority")
})

test_that("parse_karyo() errors when rules is a plain data frame", {
  df <- as.data.frame(myeloid_rules)
  class(df) <- setdiff(class(df), "karyo_rules")
  expect_error(parse_karyo("46,XX", rules = df), "karyo_rules")
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
