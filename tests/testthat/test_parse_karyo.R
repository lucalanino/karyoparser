# Helper: parse with default test options
pk <- function(...) parse_karyo(..., on_issues = "warn", verbose = FALSE)

# =============================================================================
# 1. Specific translocations (priority 100)
# =============================================================================

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

# =============================================================================
# 2. Variant band rules (priority 95)
# =============================================================================

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

# =============================================================================
# 3. Variable partner rules (priority 90)
# =============================================================================

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

# =============================================================================
# 4. Chromosome-specific rules (priority 85)
# =============================================================================

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
  # Regex idic\(X\)\(q13 is a prefix match; strip_bands() removes .1 sub-band
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
  # p-arm deletion should not match q-arm rule
  r <- pk("46,XY,del(7)(p11)")
  expect_equal(r$`del(7q)`, 0L)
})

test_that("del(12)(q22) does NOT trigger del(12p)", {
  # q-arm deletion should not match p-arm rule
  r <- pk("46,XX,del(12)(q22)")
  expect_equal(r$`del(12p)`, 0L)
})

# =============================================================================
# 5. General rules (priority 60-80)
# =============================================================================

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
  # Regex ^[?~]?t\([0-9XY] explicitly includes Y
  r <- pk("46,XY,t(Y;1)(q11;p11)")
  expect_equal(r$general_translocation, 1L)
})

# =============================================================================
# 6. Priority system
# =============================================================================

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

# =============================================================================
# 7. Ploidy
# =============================================================================

test_that("ploidy_from_count covers all thresholds", {
  expect_equal(ploidy_from_count(25), "near_haploid")
  expect_equal(ploidy_from_count(23), "near_haploid")
  expect_equal(ploidy_from_count(29), "near_haploid")
  expect_equal(ploidy_from_count(30), "low_hypodiploid")
  expect_equal(ploidy_from_count(33), "low_hypodiploid")
  expect_equal(ploidy_from_count(35), "other") # gap 34-39
  expect_equal(ploidy_from_count(40), "high_hypodiploid")
  expect_equal(ploidy_from_count(45), "high_hypodiploid")
  expect_equal(ploidy_from_count(46), "diploid")
  expect_equal(ploidy_from_count(48), "other") # gap 47-50
  expect_equal(ploidy_from_count(51), "hyperdiploid")
  expect_equal(ploidy_from_count(92), "hyperdiploid")
  expect_equal(ploidy_from_count(NA), "unknown")
})

test_that("ploidy_category for composite karyotype picks most abnormal", {
  # 46 (diploid) vs 42 (high_hypodiploid) -> should pick 42
  res <- ploidy_category("46,XX[10]/42,XX,-5,-7,-8,-9[10]")
  expect_equal(res$ploidy, "high_hypodiploid")
  expect_equal(res$chromosome_count, 42L)
})

test_that("ploidy_category filters clones with <5 metaphases", {
  # Clone with 3 metaphases should be excluded; clone with 10 used
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
  # [3~5] -> max = 5, so clone is eligible (>= 5); [3] clone is excluded.
  # 47 falls in the 47-50 "other" gap; what matters is that the 47-count clone
  # was selected (not the excluded 25-count clone), confirming range-max logic.
  res <- ploidy_category("25,X[3]/47,XX,+8[3~5]")
  expect_equal(res$ploidy, "other")
  expect_equal(res$chromosome_count, 47L)
})

test_that("ploidy_from_count: values below near_haploid floor return 'other'", {
  expect_equal(ploidy_from_count(22), "other") # below near_haploid (23-29)
})

test_that("ploidy_from_count: exact gap boundaries return 'other'", {
  # 34-39 gap (between low_hypodiploid ceiling 33 and high_hypodiploid floor 40)
  expect_equal(ploidy_from_count(34), "other")
  expect_equal(ploidy_from_count(39), "other")
  # 47-50 gap (between diploid 46 and hyperdiploid floor 51)
  expect_equal(ploidy_from_count(47), "other")
  expect_equal(ploidy_from_count(50), "other")
})

test_that("ploidy_category: all clones with <5 metaphases returns 'unknown'", {
  # has_any_brackets=TRUE but no clone meets the >=5 metaphase threshold
  res <- ploidy_category("46,XX[2]/45,XY,-7[3]")
  expect_equal(res$ploidy, "unknown")
  expect_true(is.na(res$chromosome_count))
  expect_false(res$mixed)
})

test_that("ploidy_category: range chromosome count uses round(mean())", {
  # 45~46 -> mean=45.5 -> round(45.5)=46 (R banker's rounding) -> "diploid"
  # No brackets -> all clones eligible
  res <- ploidy_category("45~46,XX")
  expect_equal(res$ploidy, "diploid")
  expect_equal(res$chromosome_count, 46L)
})

# =============================================================================
# 8. Monosomy / trisomy
# =============================================================================

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

# =============================================================================
# 9. Complex karyotype
# =============================================================================

test_that("complex karyotype threshold at 3 unique aberrations", {
  # 2 aberrations -> not complex
  r2 <- pk("46,XX,del(5)(q13),del(7)(q22)")
  expect_equal(r2$complex_karyotype, 0L)
  # 3 aberrations -> complex
  r3 <- pk("46,XX,del(5)(q13),del(7)(q22),+8")
  expect_equal(r3$complex_karyotype, 1L)
})

test_that("idem and sl excluded from complex count", {
  # Clone 1 has 2 aberrations, clone 2 adds idem + 1 new = 3 unique total (>= 3 threshold)
  r <- pk("46,XX,del(5)(q13),del(7)(q22)[10]/46,idem,+8[5]")
  expect_equal(r$complex_karyotype, 1L)
  # idem itself doesn't count as an aberration; only unique non-idem/sl tokens are counted
})

# =============================================================================
# 10. Monosomal karyotype
# =============================================================================

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
  expect_equal(r$derivative_chromosome, 1L)
})

test_that("+der() detected as derivative_chromosome", {
  r <- pk("47,XX,+der(21)")
  expect_equal(r$derivative_chromosome, 1L)
})

test_that("autosomal mono + +der() -> monosomal", {
  r <- pk("45,XX,-7,+der(21)")
  expect_equal(r$derivative_chromosome, 1L)
  expect_equal(r$monosomal_karyotype, 1L)
})

test_that("derivative_chromosome: autosomal mono + der() -> monosomal", {
  r <- pk("45,XX,-7,der(1)t(2;3)(p11;q22)")
  expect_equal(r$derivative_chromosome, 1L)
  expect_equal(r$monosomal_karyotype, 1L)
})

test_that("derivative_chromosome: co-fires with specific rule for same token", {
  # der(5)t(5;8) — t(5q) fires (translocation group) AND derivative_chromosome
  # fires (derivative group); different groups => both fire
  r <- pk("46,XX,der(5)t(5;8)(q11;q11)")
  expect_equal(r$`t(5q)`, 1L)
  expect_equal(r$derivative_chromosome, 1L)
})

test_that("derivative_chromosome: CBF-AML der co-fires but monosomal is overridden", {
  # der(8)t(8;21) — t(8;21) fires (translocation) AND derivative_chromosome fires
  # (derivative); CBF override reverts monosomal_karyotype to 0
  r <- pk("45,XX,-7,der(8)t(8;21)(q22;q22)")
  expect_equal(r$`t(8;21)(q22;q22)`, 1L)
  expect_equal(r$derivative_chromosome, 1L)
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("ider() detected as derivative_chromosome", {
  r <- pk("46,XX,ider(1)(q10)t(1;4)(p22;q11)")
  expect_equal(r$derivative_chromosome, 1L)
})

test_that("1 autosomal mono alone -> not monosomal", {
  r <- pk("45,XY,-7")
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("sex chromosome monosomy alone doesn't qualify", {
  r <- pk("45,X,-Y")
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("marker_chromosome does NOT count for monosomal", {
  r <- pk("46,XY,-7,+mar")
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("sex chromosome monosomy + structural does NOT qualify as monosomal", {
  # -X is sex chromosome monosomy; only autosomal monosomies count
  r <- pk("45,X,del(5)(q13)")
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("CBF-AML inv(16) overrides monosomal karyotype to 0", {
  # Same logic as t(8;21): CBF-AML translocations suppress monosomal flag
  r <- pk("45,XX,-7,inv(16)(p13q22)")
  expect_equal(r$`inv(16)(p13q22)`, 1L)
  expect_equal(r$monosomal_karyotype, 0L)
})

test_that("CBF-AML t(16;16) overrides monosomal karyotype to 0", {
  r <- pk("45,XX,-7,t(16;16)(p13;q22)")
  expect_equal(r$`t(16;16)(p13;q22)`, 1L)
  expect_equal(r$monosomal_karyotype, 0L)
})

# =============================================================================
# 11. Preprocessing
# =============================================================================

test_that("unicode normalization: NBSPs detected and fixed by preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo("46,\u00A0XX"))
  expect_equal(result$preprocessed, "46,XX")
  expect_equal(result$status, "fixed")
  # on_issues="warn": unicode_notation is a flagged issue — row returns NA
  r_warn <- pk("46,\u00A0XX")
  expect_true(is.na(r_warn$normal_karyotype))
  # on_issues="fix": preprocessed through .assess_karyotypes — parses correctly
  r_fix <- parse_karyo("46,\u00A0XX", on_issues = "fix", verbose = FALSE)
  expect_equal(r_fix$normal_karyotype, 1L)
})

test_that("whitespace and delimiter cleanup: normalize_iscn runs inside preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo("46 , XX , +8"))
  expect_equal(result$preprocessed, "46,XX,+8")
  # on_issues="fix" parses the normalized string
  r <- parse_karyo("46 , XX , +8", on_issues = "fix", verbose = FALSE)
  expect_equal(r$tris8, 1L)
})

test_that("psu dic normalization: preprocess_karyo normalizes before parse", {
  result <- suppressMessages(preprocess_karyo("46,XX,PSU DIC(15;22)"))
  expect_equal(result$preprocessed, "46,XX,psu dic(15;22)")
  r <- parse_karyo("46,XX,PSU DIC(15;22)", on_issues = "fix", verbose = FALSE)
  expect_equal(r$pseudodicentric, 1L)
})

test_that("idem and sl case normalization via preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo(
    "46,XX,del(5)(q13)[10]/46,IDEM,+8[5]"
  ))
  expect_equal(result$preprocessed, "46,XX,del(5)(q13)[10]/46,idem,+8[5]")
  r <- parse_karyo(
    "46,XX,del(5)(q13)[10]/46,IDEM,+8[5]",
    on_issues = "fix",
    verbose = FALSE
  )
  expect_equal(r$comma_count_aberrations, 2L)
})

test_that("duplicate commas and leading/trailing delimiters: normalize_iscn in preprocess", {
  result <- suppressMessages(preprocess_karyo(",46,,XX,"))
  expect_equal(result$preprocessed, "46,XX")
  r <- parse_karyo(",46,,XX,", on_issues = "fix", verbose = FALSE)
  expect_equal(r$normal_karyotype, 1L)
})

test_that("fullwidth plus: detected as unicode_notation, fixed by preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo("47,XY,\uFF0B8"))
  expect_equal(result$preprocessed, "47,XY,+8")
  r <- parse_karyo("47,XY,\uFF0B8", on_issues = "fix", verbose = FALSE)
  expect_equal(r$tris8, 1L)
})

test_that("cp bracket spacing normalization via parse_karyo", {
  r <- pk("46,XX[cp 20]")
  expect_true(!is.na(r$total_metaphases))
})

test_that("range count: dash normalized to tilde via parse_karyo", {
  r <- pk("46-47,XY,+8")
  expect_equal(r$normal_karyotype, 0L) # range karyotype is never normal
})

test_that("karyotype starting with 45 correctly assigns monosomy", {
  r <- pk("45,XY,-7")
  expect_equal(r$mono7, 1L)
})

# =============================================================================
# 12. Idem expansion
# =============================================================================

test_that("idem inherits stemline comma count", {
  r <- pk("46,XX,del(5)(q13),del(7)(q22)[10]/46,idem,+8[5]")
  expect_equal(r$comma_count_aberrations, 3L)
})

test_that("clone without idem uses raw count", {
  r <- pk("46,XX,del(5)(q13)")
  expect_equal(r$comma_count_aberrations, 1L)
})

test_that("idem in clone 2 inheriting 0 aberrations from clean stemline", {
  # Clone 1 has no aberrations; clone 2 inherits via idem -> 0 count
  # invalid_idem must NOT fire (idem is in clone 2, not clone 1)
  result <- suppressMessages(check_karyo("46,XX[10]/46,idem[5]"))
  expect_equal(result$invalid_idem, 0L)
  r <- pk("46,XX[10]/46,idem[5]")
  expect_false(is.na(r$ploidy_category))
  expect_equal(r$comma_count_aberrations, 0L)
  expect_equal(r$complex_karyotype, 0L)
})

# =============================================================================
# 13. check_karyo()
# =============================================================================

test_that("check_karyo: clean input returns one row per string, all zeros", {
  result <- suppressMessages(check_karyo(c("46,XX", "47,XY,+21[10]")))
  expect_equal(nrow(result), 2L)
  expect_equal(result$fixable, c(0L, 0L))
  expect_equal(result$unfixable, c(0L, 0L))
})

test_that("check_karyo: column schema matches .all_issue_types", {
  result <- suppressMessages(check_karyo("46,XX"))
  expected_cols <- c(
    "row_index",
    "karyotype",
    karyoparser:::.all_issue_types,
    "fixable",
    "unfixable"
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

test_that("check_karyo: dirty markers detected as fixable", {
  result <- suppressMessages(check_karyo(".47,XY,+21"))
  expect_equal(result$leading_dot, 1L)
  expect_equal(result$fixable, 1L)
  expect_equal(result$unfixable, 0L)
})

test_that("check_karyo: trailing narrative detected as fixable", {
  result <- suppressMessages(check_karyo("46,XX[20] .some text"))
  expect_equal(result$trailing_narrative, 1L)
  expect_equal(result$fixable, 1L)
})

test_that("check_karyo: html entities detected as fixable", {
  result <- suppressMessages(check_karyo("46,XX,t(9;22)(q34;q11) &lt;AML&gt;"))
  expect_equal(result$html_entities, 1L)
  expect_equal(result$fixable, 1L)
})

test_that("check_karyo: one row per input; clean/dirty rows correctly flagged", {
  x <- c("46,XX", ".47,XY,+21", "46,XX[5] .note", "46,XY")
  result <- suppressMessages(check_karyo(x))
  expect_equal(nrow(result), 4L)
  expect_equal(result$fixable, c(0L, 1L, 1L, 0L))
})

test_that("check_karyo: multiple issues on same row all flagged", {
  result <- suppressMessages(check_karyo(".46,XX[10] .note"))
  expect_equal(result$leading_dot, 1L)
  expect_equal(result$trailing_narrative, 1L)
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
  expect_false(is.na(r$ploidy_category))
  expect_equal(r$tris1, 1L)
})

test_that("check_karyo: no_sex_complement detected when sex chromosome token absent", {
  result <- suppressMessages(check_karyo("46,+8"))
  expect_equal(result$no_sex_complement, 1L)
  expect_equal(result$unfixable, 1L)
})

# =============================================================================
# 14. on_issues parameter
# =============================================================================

test_that("on_issues='warn': issue rows get NA output", {
  dirty <- c("46,XX", ".47,XY,+21")
  result <- pk(dirty)
  expect_equal(nrow(result), 2L)
  # Row 1 (clean) should have a ploidy_category
  expect_false(is.na(result$ploidy_category[1]))
  # Row 2 (dirty) should be NA
  expect_true(is.na(result$ploidy_category[2]))
})

test_that("on_issues='warn': emits message for dirty input", {
  dirty <- c("46,XX", ".47,XY,+21")
  expect_message(
    parse_karyo(dirty, on_issues = "warn", verbose = TRUE),
    regexp = "issues"
  )
})

test_that("on_issues='stop': throws error for dirty input", {
  dirty <- c("46,XX", ".47,XY,+21")
  expect_error(
    parse_karyo(dirty, on_issues = "stop", verbose = FALSE),
    regexp = "have issues"
  )
})

test_that("on_issues='stop': throws error for invalid input", {
  expect_error(
    parse_karyo(c("46,XX", NA), on_issues = "stop", verbose = FALSE),
    regexp = "have issues"
  )
})

test_that("on_issues='warn': invalid rows get NA", {
  r <- pk(c("46,XX", NA))
  expect_equal(r$normal_karyotype[1], 1L)
  expect_true(is.na(r$normal_karyotype[2]))
})

test_that("on_issues: clean input produces no guard message", {
  clean <- c("46,XX", "47,XY,+21[10]")
  expect_no_message(
    parse_karyo(clean, on_issues = "warn", verbose = FALSE)
  )
})

# =============================================================================
# 15. fixable_error and unfixable_error columns
# =============================================================================

test_that("error columns: clean row has both = 0", {
  r <- pk("46,XX")
  expect_equal(r$fixable_error, 0L)
  expect_equal(r$unfixable_error, 0L)
})

test_that("error columns: dirty row under on_issues='warn' has fixable_error=1, unfixable_error=0", {
  r <- pk(c("46,XX", ".47,XY,+21"))
  expect_equal(r$fixable_error[1], 0L)
  expect_equal(r$fixable_error[2], 1L)
  expect_equal(r$unfixable_error[2], 0L)
})

test_that("error columns: structural error row has unfixable_error=1, fixable_error=0", {
  r <- pk(c("46,XX", NA))
  expect_equal(r$unfixable_error[1], 0L)
  expect_equal(r$unfixable_error[2], 1L)
  expect_equal(r$fixable_error[2], 0L)
})

test_that("error columns: successfully fixed row has both = 0", {
  r <- parse_karyo(c("46,XX", ".47,XY,+21"), on_issues = "fix", verbose = FALSE)
  expect_equal(r$fixable_error[2], 0L)
  expect_equal(r$unfixable_error[2], 0L)
})

test_that("error columns: present in empty result", {
  r <- pk(character(0))
  expect_true("fixable_error" %in% names(r))
  expect_true("unfixable_error" %in% names(r))
})

test_that("error columns: survive dplyr::filter()", {
  r <- pk(c("46,XX", NA, "47,XY,+8"))
  filtered <- dplyr::filter(r, !is.na(normal_karyotype))
  expect_equal(nrow(filtered), 2L)
  expect_true("fixable_error" %in% names(filtered))
  expect_true("unfixable_error" %in% names(filtered))
})

# =============================================================================
# 16. ID column
# =============================================================================

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
  r <- pk(df)
  expect_equal(r$sample_id, "A")
  expect_true(is.na(r$ploidy_category))
  expect_equal(names(r)[1], "sample_id")
})

test_that("explicit karyotype_column parameter works for non-standard column name", {
  df <- data.frame(my_col = c("46,XX", "46,XY"), stringsAsFactors = FALSE)
  r <- pk(df, karyotype_column = "my_col")
  expect_equal(nrow(r), 2L)
  expect_equal(r$normal_karyotype, c(1L, 1L))
})

# =============================================================================
# 17. Normal karyotype
# =============================================================================

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

# =============================================================================
# 18. Integration tests
# =============================================================================

test_that("full pipeline: normal karyotype", {
  r <- pk("46,XX")
  expect_equal(r$ploidy_category, "diploid")
  expect_equal(r$normal_karyotype, 1L)
  expect_equal(r$complex_karyotype, 0L)
  expect_equal(r$monosomal_karyotype, 0L)
  expect_equal(r$chromosome_count, 46L)
})

test_that("full pipeline: APL t(15;17)", {
  r <- pk("46,XX,t(15;17)(q24;q21)")
  expect_equal(r$`t(15;17)(q24;q21)`, 1L)
  expect_equal(r$ploidy_category, "diploid")
  expect_equal(r$normal_karyotype, 0L)
  expect_equal(r$complex_karyotype, 0L)
})

test_that("full pipeline: complex monosomal karyotype", {
  r <- pk("43,XY,-5,-7,del(17)(p13),+8")
  expect_equal(r$mono5, 1L)
  expect_equal(r$mono7, 1L)
  expect_equal(r$tris8, 1L)
  expect_equal(r$`del(17p)`, 1L)
  expect_equal(r$complex_karyotype, 1L)
  expect_equal(r$monosomal_karyotype, 1L)
})

test_that("full pipeline: composite karyotype with metaphases", {
  r <- pk("46,XX,del(5)(q13)[10]/47,XX,del(5)(q13),+8[5]")
  expect_true(!is.na(r$total_metaphases))
  expect_equal(r$total_metaphases, 15L)
  expect_equal(r$`del(5q)`, 1L)
  expect_equal(r$tris8, 1L)
})

test_that("full pipeline: hyperdiploid karyotype", {
  r <- pk("51,XY,+8,+11,+13,+19,+21")
  expect_equal(r$ploidy_category, "hyperdiploid")
  expect_equal(r$chromosome_count, 51L)
  expect_equal(r$tris8, 1L)
  expect_equal(r$tris21, 1L)
})

# =============================================================================
# 19. Edge cases
# =============================================================================

test_that("empty input returns empty tibble", {
  r <- pk(character(0))
  expect_equal(nrow(r), 0)
  expect_true("original_karyotype" %in% names(r))
})

test_that("single karyotype input works", {
  r <- pk("46,XX")
  expect_equal(nrow(r), 1)
})

test_that("version attribute is set", {
  r <- pk("46,XX")
  expect_equal(attr(r, "karyoparser_version"), "0.7.1")
})

test_that(".return='data.frame' returns data.frame", {
  r <- parse_karyo(
    "46,XX",
    .return = "data.frame",
    on_issues = "warn",
    verbose = FALSE
  )
  expect_true(is.data.frame(r))
  expect_false(tibble::is_tibble(r))
})

test_that("comma_count_aberrations for simple karyotype", {
  # 46,XX -> sex complement only, 0 aberrations
  r <- pk("46,XX")
  expect_equal(r$comma_count_aberrations, 0L)
})

test_that("total_metaphases is NA when no brackets", {
  r <- pk("46,XX")
  expect_true(is.na(r$total_metaphases))
})

test_that("mixed_ploidy flag in output", {
  r <- pk("45,XX,-7[10]/51,XX,+8,+11,+13,+19,+21[8]")
  expect_equal(r$mixed_ploidy, 1L)
})

# =============================================================================
# 20. Deduplication
# =============================================================================

test_that("duplicate karyotypes produce same results as unique input", {
  input <- c("46,XX,del(5)(q13)", "47,XY,+8", "46,XX,del(5)(q13)")
  r <- pk(input)
  expect_equal(nrow(r), 3)
  # Rows 1 and 3 should be identical (same karyotype)
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
  expect_equal(r$`del(7q)`[1], r$`del(7q)`[3])
})

test_that("dedup handles mix of valid and invalid karyotypes", {
  input <- c("46,XX", NA, "46,XX", "47,XY,+8", NA)
  r <- pk(input)
  expect_equal(nrow(r), 5)
  expect_equal(r$normal_karyotype[1], 1L)
  expect_equal(r$normal_karyotype[3], 1L)
  expect_true(is.na(r$normal_karyotype[2]))
  expect_true(is.na(r$normal_karyotype[5]))
})

# =============================================================================
# 21. Multi-group rule firing
# =============================================================================

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

# =============================================================================
# 22. preprocess_karyo()
# =============================================================================

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

test_that("preprocess_karyo: decodes &amp; entity", {
  expect_equal(
    suppressMessages(preprocess_karyo("46,XX &amp; notes"))$preprocessed,
    "46,XX & notes"
  )
})

test_that("preprocess_karyo: .// and // prefix are unfixable; normalize_iscn strips the prefix from output", {
  result <- suppressMessages(preprocess_karyo(c(".//46,XX", "..//46,XX")))
  expect_equal(result$preprocessed, c("46,XX", "46,XX"))
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

# =============================================================================
# 23. on_issues ("fix"/"warn"/"stop"); .dirty_patterns
# =============================================================================

# 23a. on_issues = "fix" ---------------------------------------------------

test_that("on_issues='fix': dirty row is cleaned and parsed (not NA)", {
  r <- parse_karyo(
    c("46,XX", ".47,XY,+21"),
    on_issues = "fix",
    verbose = FALSE
  )
  expect_equal(nrow(r), 2L)
  expect_equal(r$tris21[2], 1L)
  expect_false(is.na(r$ploidy_category[2]))
})

test_that("on_issues='fix': structural errors still produce NA", {
  r <- parse_karyo(
    c("46,XX", "not_a_karyotype"),
    on_issues = "fix",
    verbose = FALSE
  )
  expect_true(is.na(r$ploidy_category[2]))
})

test_that("on_issues='fix': verbose message reports parsed/fixed counts", {
  expect_message(
    parse_karyo(c("46,XX", ".47,XY,+21"), on_issues = "fix", verbose = TRUE),
    "fixed"
  )
})

test_that("on_issues='fix': successfully cleaned row has both error columns = 0", {
  r <- parse_karyo(
    c("46,XX", ".47,XY,+21"),
    on_issues = "fix",
    verbose = FALSE
  )
  expect_equal(r$fixable_error[2], 0L)
  expect_equal(r$unfixable_error[2], 0L)
})

test_that("on_issues='fix': multiple dirty rows all cleaned", {
  r <- parse_karyo(
    c(".46,XX", ".47,XY,+21"),
    on_issues = "fix",
    verbose = FALSE
  )
  expect_equal(r$normal_karyotype[1], 1L)
  expect_equal(r$tris21[2], 1L)
})

# 23b. on_issues = "warn" ------------------------------------------------------

test_that("on_issues='warn': dirty row returns NA", {
  r <- parse_karyo(
    c("46,XX", ".47,XY,+21"),
    on_issues = "warn",
    verbose = FALSE
  )
  expect_true(is.na(r$ploidy_category[2]))
})

test_that("on_issues='warn': clean row in same batch unaffected", {
  r <- parse_karyo(
    c("46,XX", ".47,XY,+21"),
    on_issues = "warn",
    verbose = FALSE
  )
  expect_equal(r$normal_karyotype[1], 1L)
})

# 23c. .dirty_patterns consistency -------------------------------------------

test_that(".dirty_patterns contains all expected keys", {
  expect_setequal(
    names(karyoparser:::.dirty_patterns),
    c(
      "unicode_notation",
      "embedded_newline",
      "html_entities",
      "leading_dot",
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
  expect_equal(
    suppressMessages(preprocess_karyo("&lt;46&gt;,XX"))$preprocessed,
    "<46>,XX"
  )
  expect_equal(
    suppressMessages(preprocess_karyo(".46,XY,+21[10] .Note"))$preprocessed,
    "46,XY,+21[10]"
  )
  # .// prefix is unfixable; normalize_iscn strips it from the preprocessed output
  expect_equal(
    suppressMessages(preprocess_karyo(".//46,XX"))$preprocessed,
    "46,XX"
  )
})

# 23d. missing_sex_comma ------------------------------------------------------

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
  res <- parse_karyo(k, on_issues = "fix", verbose = FALSE)
  expect_equal(res$ploidy_category, "diploid")
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

# =============================================================================
# 24. trailing_narrative v2, midstring_linewrap, chimeric_separator, updated_iscn
# =============================================================================

# 24a. trailing_narrative anchor rule ----------------------------------------

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

# 24b. ] ./ separator artifact -----------------------------------------------

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

# 24c. midstring_linewrap ----------------------------------------------------

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
  # ' .Capital' is trailing narrative, not mid-string linewrap
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX[20] .Female karyotype"
    ))$preprocessed,
    "46,XX[20]"
  )
})

# 24d. chimeric_separator detection ------------------------------------------

test_that("check_karyo detects chimeric_separator as fixable", {
  result <- suppressMessages(check_karyo("46,XX[15]//46,XY[5]"))
  expect_equal(result$chimeric_separator, 1L)
  expect_equal(result$fixable, 1L)
})

test_that("check_karyo: clean karyotype does not trigger chimeric_separator", {
  result <- suppressMessages(check_karyo("46,XX[20]"))
  expect_equal(result$chimeric_separator, 0L)
})

# 24e. updated_iscn detection ------------------------------------------------

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
  r <- pk("46,XX Updated ISCN new")
  expect_true(is.na(r$ploidy_category))
})

test_that("updated_iscn rows stay NA even with on_issues='fix'", {
  r <- parse_karyo("46,XX Updated ISCN new", on_issues = "fix", verbose = FALSE)
  expect_true(is.na(r$ploidy_category))
})

# 24f. chimeric: on_issues = "warn" ------------------------------------------

test_that("on_issues='warn': chimeric row returns NA", {
  r <- parse_karyo("46,XX[15]//46,XY[5]", on_issues = "warn", verbose = FALSE)
  expect_true(is.na(r$ploidy_category))
})

test_that("on_issues='warn': chimeric row has fixable_error=1", {
  r <- parse_karyo("46,XX[15]//46,XY[5]", on_issues = "warn", verbose = FALSE)
  expect_equal(r$fixable_error, 1L)
  expect_equal(r$unfixable_error, 0L)
})

# 24g. chimeric: on_issues = "fix" (truncate) ---------------------------------

test_that("on_issues='fix': chimeric row truncated and parsed", {
  r <- parse_karyo(
    "46,XX,t(9;22)(q34;q11)[15]//46,XX[5]",
    on_issues = "fix",
    verbose = FALSE
  )
  expect_false(is.na(r$ploidy_category))
  expect_equal(r$`t(9;22)(q34;q11)`, 1L)
})

test_that("on_issues='fix': clean row in same batch unaffected by chimeric truncation", {
  r <- parse_karyo(
    c("46,XX[20]", "46,XX,t(9;22)(q34;q11)[15]//46,XX[5]"),
    on_issues = "fix",
    verbose = FALSE
  )
  expect_equal(r$normal_karyotype[1], 1L)
  expect_false(is.na(r$ploidy_category[2]))
})

test_that("on_issues='fix': verbose message reports fixed count for chimeric row", {
  expect_message(
    parse_karyo(
      "46,XX,t(9;22)(q34;q11)[15]//46,XX[5]",
      on_issues = "fix",
      verbose = TRUE
    ),
    "fixed"
  )
})

test_that("on_issues='fix': chimeric row with residual structural issue after truncation becomes NA", {
  r <- parse_karyo(
    "not_valid//46,XX[5]",
    on_issues = "fix",
    verbose = FALSE
  )
  expect_true(is.na(r$ploidy_category))
})

# 24h. chimeric: on_issues = "stop" ------------------------------------------

test_that("on_issues='stop': raises error on chimeric row", {
  expect_error(
    parse_karyo("46,XX[15]//46,XY[5]", on_issues = "stop", verbose = FALSE),
    "chimeric"
  )
})

# 24i. zero_host_chimera -------------------------------------------------------

test_that("check_karyo detects zero_host_chimera as unfixable", {
  result <- suppressMessages(check_karyo(".//46,XX[10]"))
  expect_equal(result$zero_host_chimera, 1L)
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

test_that("parse_karyo: zero_host_chimera returns NA", {
  r <- pk(".//46,XX[10]")
  expect_true(is.na(r$ploidy_category))
})

test_that("parse_karyo: zero_host_chimera has unfixable_error=1, fixable_error=0", {
  r <- pk(".//46,XX[10]")
  expect_equal(r$unfixable_error, 1L)
  expect_equal(r$fixable_error, 0L)
})

test_that("parse_karyo: zero_host_chimera NA even with on_issues='fix'", {
  r <- parse_karyo(".//46,XX[10]", on_issues = "fix", verbose = FALSE)
  expect_true(is.na(r$ploidy_category))
})

test_that("parse_karyo: zero_host_chimera with other dirty patterns still returns NA", {
  # Trailing narrative present alongside .// — other dirty fix runs but row is still NA
  r <- parse_karyo(
    ".//46,XX[10] .Female karyotype",
    on_issues = "fix",
    verbose = FALSE
  )
  expect_true(is.na(r$ploidy_category))
})

test_that("parse_karyo: row with both fixable and unfixable issues has both error columns = 1", {
  # zero_host_chimera (unfixable) + trailing_narrative (fixable) co-occur
  r <- pk(".//46,XX[10] .Female karyotype")
  expect_equal(r$fixable_error, 1L)
  expect_equal(r$unfixable_error, 1L)
})

test_that("check_karyo: always prints checking count and summary", {
  expect_message(
    check_karyo(c("46,XX", ".47,XY,+21", NA)),
    "Checking 3 karyotype"
  )
  expect_message(
    check_karyo(c("46,XX", ".47,XY,+21", NA)),
    "Fixable"
  )
})

test_that("check_karyo: prints 'All clean.' when no issues", {
  expect_message(
    check_karyo("46,XX[20]"),
    "All clean"
  )
})

test_that("check_karyo: verbose=TRUE adds per-type breakdown", {
  expect_message(
    check_karyo(c("46,XX", ".47,XY,+21", NA), verbose = TRUE),
    "breakdown"
  )
})

test_that("check_karyo: data frame input errors with helpful message", {
  df <- data.frame(karyotype = "46,XX", stringsAsFactors = FALSE)
  expect_error(check_karyo(df), "data frame")
})

test_that("preprocess_karyo: data frame input errors with helpful message", {
  df <- data.frame(karyotype = "46,XX", stringsAsFactors = FALSE)
  expect_error(preprocess_karyo(df), "data frame")
})

test_that("parse_karyo: clean row in same batch unaffected by zero_host_chimera row", {
  r <- pk(c("46,XX[20]", ".//46,XY[10]"))
  expect_equal(r$normal_karyotype[1], 1L)
  expect_true(is.na(r$ploidy_category[2]))
})

# =============================================================================
# 25. fish_notation, mar_space, midstring_linewrap +
# =============================================================================

# 25a. fish_notation -----------------------------------------------------------

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

# 25b. mar_space --------------------------------------------------------------

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

# 25c. midstring_linewrap with + ----------------------------------------------

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

# =============================================================================
# 26. B2/B4/B5/chimeric_karyotype (v0.6.0)
# =============================================================================

# 26a. B2: fix handles both dot and non-dot linewrap variants -----------------
# Note: detect still requires a dot (to avoid false positives on "46 , XX , +8").
# The fix runs unconditionally so bare-space wraps are still cleaned.

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
  # Not detected as midstring_linewrap (no dot), but normalize_iscn in preprocess tidies it
  result <- suppressMessages(
    preprocess_karyo("46,XY,del(5)(q13), t(9;22)(q34;q11.2)[10]/46,XY[5]")
  )
  expect_equal(
    result$preprocessed,
    "46,XY,del(5)(q13),t(9;22)(q34;q11.2)[10]/46,XY[5]"
  )
})

# 26b. B4: trailing narrative without metaphase bracket -----------------------

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
  # Rule 1 handles this; rule 3 must not double-strip
  expect_equal(
    suppressMessages(preprocess_karyo(
      "46,XX,t(9;22)(q34;q11.2)[20] Abnormal female karyotype"
    ))$preprocessed,
    "46,XX,t(9;22)(q34;q11.2)[20]"
  )
})

# 26c. B5: space before '[' normalized ----------------------------------------

test_that("normalize_iscn: collapses space before opening bracket", {
  result <- suppressMessages(pk("46,XX [20]"))
  expect_equal(result$normal_karyotype, 1L)
})

test_that("preprocess_karyo: space before '[' cleaned via normalize_iscn in preprocess_karyo", {
  result <- suppressMessages(preprocess_karyo("46,XX [20]"))
  expect_equal(result$preprocessed, "46,XX[20]")
})

# 26d. chimeric_karyotype column -----------------------------------------------

test_that("parse_karyo: chimeric_karyotype column exists in output", {
  result <- suppressMessages(pk("46,XX"))
  expect_true("chimeric_karyotype" %in% names(result))
})

test_that("parse_karyo: clean row has chimeric_karyotype=0", {
  result <- suppressMessages(pk("46,XX"))
  expect_equal(result$chimeric_karyotype, 0L)
})

test_that("parse_karyo on_issues='fix': chimeric row has chimeric_karyotype=1", {
  result <- suppressMessages(
    parse_karyo(
      "46,XX,t(9;22)(q34;q11.2)[3]//46,XY[12]",
      on_issues = "fix",
      verbose = FALSE
    )
  )
  expect_equal(result$chimeric_karyotype, 1L)
})

test_that("parse_karyo on_issues='warn': chimeric row has chimeric_karyotype=1", {
  result <- suppressMessages(
    parse_karyo(
      "46,XX,t(9;22)(q34;q11.2)[3]//46,XY[12]",
      on_issues = "warn",
      verbose = FALSE
    )
  )
  expect_equal(result$chimeric_karyotype, 1L)
})

test_that("parse_karyo: zero_host_chimera row has chimeric_karyotype=1", {
  result <- suppressMessages(
    parse_karyo(".//46,XX[20]", on_issues = "fix", verbose = FALSE)
  )
  expect_equal(result$chimeric_karyotype, 1L)
})

test_that("parse_karyo: zero_host_chimera has both unfixable_error=1 and chimeric_karyotype=1", {
  result <- suppressMessages(
    parse_karyo(".//46,XX[20]", on_issues = "fix", verbose = FALSE)
  )
  expect_equal(result$unfixable_error, 1L)
  expect_equal(result$chimeric_karyotype, 1L)
})

test_that("parse_karyo: non-chimeric row in chimeric batch has chimeric_karyotype=0", {
  result <- suppressMessages(
    parse_karyo(
      c("46,XX,t(9;22)(q34;q11.2)[3]//46,XY[12]", "46,XY"),
      on_issues = "fix",
      verbose = FALSE
    )
  )
  expect_equal(result$chimeric_karyotype, c(1L, 0L))
})

# =============================================================================
# 27. API surface: rules_table(), normalized_karyotype
# =============================================================================

test_that("rules_table() returns a tibble with expected columns", {
  rt <- rules_table()
  expect_true(tibble::is_tibble(rt))
  expect_true(nrow(rt) > 0L)
  expect_true(all(
    c(
      "flag_name",
      "regex",
      "category",
      "priority",
      "counts_for_monosomal",
      "competition_group"
    ) %in%
      names(rt)
  ))
})

test_that("rules_table() priorities are positive integers", {
  rt <- rules_table()
  expect_true(is.integer(rt$priority) || is.numeric(rt$priority))
  expect_true(all(rt$priority >= 1L))
})

test_that("normalized_karyotype under on_issues='warn' equals the raw input", {
  r <- parse_karyo(c("46,XX", "46,XY,+8"), on_issues = "warn", verbose = FALSE)
  expect_equal(r$normalized_karyotype[1], "46,XX")
  expect_equal(r$normalized_karyotype[2], "46,XY,+8")
})

test_that("normalized_karyotype under on_issues='fix' equals the cleaned string", {
  r <- parse_karyo(".46,XX", on_issues = "fix", verbose = FALSE)
  expect_equal(r$normalized_karyotype, "46,XX")
})

test_that("normalized_karyotype column is present in output", {
  r <- pk("46,XX")
  expect_true("normalized_karyotype" %in% names(r))
})
