#' Assign Cytogenetic Risk Categories
#'
#' Adds a cytogenetic risk category to each row of a [parse_karyo()] result.
#' Interpretation only -- no parsing or cleaning happens here; the function
#' reads columns [parse_karyo()] already produced.
#'
#' **This is the cytogenetic component only, never a complete risk score.**
#' IPSS-R stratifies on five components -- cytogenetics, marrow blast
#' percentage, haemoglobin, platelet count, and absolute neutrophil count --
#' of which this computes the first. ELN 2022 is a *genetic* classification:
#' assigning it properly needs *NPM1*, *FLT3*-ITD, *CEBPA* bZIP, *TP53*,
#' *ASXL1*, *RUNX1* and others, none of which a karyotype carries. A normal
#' karyotype is intermediate cytogenetically but favorable with mutated
#' *NPM1* and no *FLT3*-ITD, or adverse with mutated *TP53*. Treat both
#' columns as one input to a risk assessment, not as the assessment.
#'
#' Apply one scheme per call; chain for both:
#' `assign_risk(x, "ipssr") |> assign_risk("eln2022")`.
#'
#' @section Counting aberrations:
#' "Isolated", "double", and "complex" are resolved with
#' `distinct_aberrations` -- distinct normalized aberration tokens pooled
#' across all clones, which is the same count `complex_karyotype` thresholds
#' at 3. Rows the parser could not score (`distinct_aberrations` is `NA`)
#' return `NA`.
#'
#' @section IPSS-R categories:
#' Following the IPSS-R cytogenetic risk groups (Schanz et al. 2012):
#'
#' | Category | Score | Definition |
#' |---|---|---|
#' | Very Good | 0 | isolated `-Y`, isolated `del(11q)` |
#' | Good | 1 | normal, isolated `del(5q)`/`del(12p)`/`del(20q)`, double including `del(5q)` |
#' | Intermediate | 2 | `del(7q)`, `+8`, `+19`, `i(17q)`, any other single or double |
#' | Poor | 3 | `-7`, `inv(3)`/`t(3q)`/`del(3q)`, double including `-7`/`del(7q)`, 3 abnormalities |
#' | Very Poor | 4 | more than 3 abnormalities |
#'
#' Two resolutions are applied where the published table is silent, both
#' deliberate and worth confirming against the source before relying on them:
#'
#' - A double containing **both** `del(5q)` and `-7`/`del(7q)` satisfies the
#'   Good row and the Poor row at once. Poor wins here, on the grounds that
#'   the more adverse lesion should not be masked by its partner.
#' - `t(3q)` is read as 3q26 rearrangement (`t_3q26_v`), not any 3q
#'   breakpoint, since the IPSS-R row reflects *MECOM*/EVI1 biology.
#'   `del(3q)` is matched at any 3q breakpoint.
#'
#' A karyotype with no scoreable aberration token that is nonetheless not
#' normal -- `45,X`, where the loss rides on the sex complement and `-X`
#' cannot be distinguished from `-Y` -- returns `NA` rather than guessing.
#'
#' @section ELN 2022 categories:
#' The cytogenetic rows of the ELN 2022 risk classification, evaluated in
#' strict precedence order:
#'
#' 1. **Favorable** -- `t(8;21)(q22;q22)`/*RUNX1::RUNX1T1*,
#'    `inv(16)(p13q22)` or `t(16;16)(p13;q22)`/*CBFB::MYH11*. These are
#'    class-defining, so they win over every adverse criterion below,
#'    mirroring the CBF-AML override already applied to
#'    `monosomal_karyotype`.
#' 2. **Intermediate, by carve-out** -- `t(9;11)(p21;q23)`/*MLLT3::KMT2A*.
#'    Checked before the adverse rows because a `t(9;11)` also satisfies the
#'    adverse `t(v;11q23)`/*KMT2A*-rearranged row; ELN gives `t(9;11)`
#'    precedence, so it must be resolved first or every `t(9;11)` would be
#'    miscalled adverse.
#' 3. **Adverse** -- `t(6;9)`/*DEK::NUP214*; `t(v;11q23)`/*KMT2A*-rearranged;
#'    `t(9;22)`/*BCR::ABL1*; `t(8;16)`/*KAT6A::CREBBP*; `inv(3)(q21q26)` or
#'    `t(3;3)(q21;q26)`; `t(3q26;v)`/*MECOM*-rearranged; `-5` or `del(5q)`;
#'    `-7`; `-17`/abn(17p); complex karyotype; monosomal karyotype.
#' 4. **Intermediate** -- anything else, including a normal karyotype.
#'
#' `t(6;9)` and `t(9;22)` are matched at both canonical and non-canonical
#' breakpoints, since ELN names the fusion rather than the breakpoint. This
#' matters: ELN writes `t(6;9)(p23.3;q34.1)`, which is a *non*-canonical
#' spelling for [myeloid_rules] and would otherwise be missed. `inv(3)` and
#' `t(3;3)` are matched only at the canonical `q21q26`, because the
#' non-canonical variants are not the *GATA2*/*MECOM* lesion.
#'
#' **Complex** uses the ELN definition, not the scheme-neutral
#' `complex_karyotype` column: three or more distinct aberrations, excluding
#' hyperdiploid karyotypes with three or more trisomies (or polysomies) and
#' no structural abnormality. The "in the absence of a class-defining
#' recurring genetic abnormality" clause falls out of the precedence order --
#' favorable lesions and `t(9;11)` are resolved before complex is considered.
#' `complex_karyotype` itself is deliberately left alone, since IPSS-R's
#' complex rows carry no such exclusion.
#'
#' Two judgment calls, both worth confirming against the source:
#'
#' - **abn(17p)** is read as `-17`, `del(17p)`, `add(17p)` or `i(17q)`.
#'   Including `i(17q)` reflects that it entails 17p loss.
#' - **Partial 17p loss implied by an unbalanced derivative is *not*
#'   counted**, because `unbal_partial_loss_*` is documented as separate from
#'   `del()`/`mono*` and folding it in here would breach that boundary. If
#'   you want those rows treated as abn(17p), combine
#'   `unbal_partial_loss_17p` yourself.
#'
#' ELN 2022 has no numeric score, so `scheme = "eln2022"` adds only a
#' category column.
#'
#' @param parsed A tibble returned by [parse_karyo()]. Must carry the default
#'   [myeloid_rules] flags; if `columns` was used to narrow the output, the
#'   classes the scheme needs must still be present.
#' @param scheme Character string naming the risk scheme: `"ipssr"` (default,
#'   IPSS-R cytogenetic risk group) or `"eln2022"` (ELN 2022 cytogenetic risk
#'   category). One per call.
#'
#' @return `parsed` with the scheme's columns appended.
#'
#'   For `"ipssr"`:
#'   - `ipssr_cyto_risk`: ordered factor, `Very Good` < `Good` <
#'     `Intermediate` < `Poor` < `Very Poor`; `NA` for rows that cannot be
#'     scored
#'   - `ipssr_cyto_score`: integer 0-4 matching the category, `NA` alongside
#'     an `NA` category
#'
#'   For `"eln2022"`:
#'   - `eln2022_cyto_risk`: ordered factor, `Favorable` < `Intermediate` <
#'     `Adverse`; `NA` for rows the parser could not score
#'
#' @examples
#' parsed <- parse_karyo(c("46,XX", "46,XX,del(5)(q13q33)", "45,XY,-7"))
#' risk <- assign_risk(parsed, scheme = "ipssr")
#' risk[, c("original_karyotype", "ipssr_cyto_risk", "ipssr_cyto_score")]
#'
#' # Both schemes, chained
#' both <- assign_risk(parsed, "ipssr") |> assign_risk("eln2022")
#' both[, c("original_karyotype", "ipssr_cyto_risk", "eln2022_cyto_risk")]
#'
#' @seealso [parse_karyo()], `vignette("risk-stratification")`
#' @export
assign_risk <- function(parsed, scheme = c("ipssr", "eln2022")) {
  scheme <- match.arg(scheme)

  if (!is.data.frame(parsed)) {
    stop(
      "`parsed` must be a data frame returned by parse_karyo().",
      call. = FALSE
    )
  }

  required <- switch(
    scheme,
    ipssr = .ipssr_required_columns,
    eln2022 = .eln2022_required_columns()
  )
  missing_cols <- setdiff(required, names(parsed))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        paste(
          "parse_karyo() output is missing column(s) required by scheme",
          "\"%s\": %s.",
          "\nRisk assignment reads the default `myeloid_rules` flags plus the",
          "aneuploidy, classification and summary columns.",
          "\nRe-run parse_karyo() with rules = myeloid_rules, and if you",
          "narrowed `columns`, keep the classes those flags belong to."
        ),
        scheme,
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  switch(
    scheme,
    ipssr = .assign_ipssr(parsed),
    eln2022 = .assign_eln2022(parsed)
  )
}

# Flags are 0/1 on parsed rows and NA on unscoreable ones; the NA rows are
# caught by `is.na(distinct_aberrations)` in each scheme, so flattening NA to
# 0 here is safe and keeps the row arithmetic total.
.risk_flag_matrix <- function(parsed, cols) {
  m <- matrix(0L, nrow = nrow(parsed), ncol = length(cols))
  for (i in seq_along(cols)) {
    m[, i] <- tidyr::replace_na(as.integer(parsed[[cols[i]]]), 0L)
  }
  m
}

.risk_count <- function(parsed, cols) {
  as.integer(rowSums(.risk_flag_matrix(parsed, cols)))
}

.risk_any <- function(parsed, cols) {
  .risk_count(parsed, cols) > 0L
}

.risk_chroms <- c(as.character(1:22), "X", "Y")

# IPSS-R ----

.ipssr_levels <- c("Very Good", "Good", "Intermediate", "Poor", "Very Poor")

.ipssr_scores <- c(
  "Very Good" = 0L,
  "Good" = 1L,
  "Intermediate" = 2L,
  "Poor" = 3L,
  "Very Poor" = 4L
)

.ipssr_required_columns <- c(
  "distinct_aberrations",
  "normal_karyotype",
  "monoY",
  "mono7",
  "del_11q",
  "del_5q",
  "del_12p",
  "del_20q",
  "del_7q",
  "inv_3_q21q26",
  "inv_3_other",
  "t_3q26_v",
  "del_3q"
)

.assign_ipssr <- function(parsed) {
  one <- function(col) .risk_any(parsed, col)

  n <- as.integer(parsed$distinct_aberrations)
  isolated <- !is.na(n) & n == 1L
  double <- !is.na(n) & n == 2L

  abn_3q <- .risk_any(
    parsed,
    c("inv_3_q21q26", "inv_3_other", "t_3q26_v", "del_3q")
  )
  mono_7 <- one("mono7")
  del_7q <- one("del_7q")
  del_5q <- one("del_5q")

  category <- dplyr::case_when(
    is.na(n) ~ NA_character_,
    n == 0L & one("normal_karyotype") ~ "Good",
    # No scoreable token but not normal (e.g. "45,X") -- see the IPSS-R
    # section of the docs.
    n == 0L ~ NA_character_,
    n > 3L ~ "Very Poor",
    n == 3L ~ "Poor",
    double & (mono_7 | del_7q) ~ "Poor",
    double & del_5q ~ "Good",
    double ~ "Intermediate",
    isolated & mono_7 ~ "Poor",
    isolated & abn_3q ~ "Poor",
    isolated & one("monoY") ~ "Very Good",
    isolated & one("del_11q") ~ "Very Good",
    isolated & del_5q ~ "Good",
    isolated & one("del_12p") ~ "Good",
    isolated & one("del_20q") ~ "Good",
    # Everything else -- including isolated del(7q), +8, +19 and i(17q), which
    # the table lists explicitly but which land on the same category as the
    # "any other single or double" catch-all.
    TRUE ~ "Intermediate"
  )

  parsed$ipssr_cyto_risk <- factor(
    category,
    levels = .ipssr_levels,
    ordered = TRUE
  )
  parsed$ipssr_cyto_score <- unname(.ipssr_scores[category])
  parsed
}

# ELN 2022 ----

.eln2022_levels <- c("Favorable", "Intermediate", "Adverse")

.eln2022_favorable <- c(
  "t_8_21_q22_q22",
  "inv_16_p13q22",
  "t_16_16_p13_q22"
)

# Both canonical and non-canonical breakpoints for t(6;9) and t(9;22): ELN
# names the fusion, not the breakpoint, and ELN's own t(6;9)(p23.3;q34.1)
# spelling is non-canonical for myeloid_rules. inv(3)/t(3;3) stay canonical,
# since their non-canonical variants are not the GATA2/MECOM lesion.
.eln2022_adverse <- c(
  "t_6_9_p22_q34",
  "t_6_9_other",
  "t_v_11q23",
  "t_9_22_q34_q11",
  "t_9_22_other",
  "t_8_16_p11_p13",
  "inv_3_q21q26",
  "t_3_3_q21_q26",
  "t_3q26_v",
  "mono5",
  "del_5q",
  "mono7",
  # abn(17p): -17, del(17p), add(17p), i(17q). Der-implied 17p loss is
  # deliberately excluded -- see the ELN section of the docs.
  "mono17",
  "del_17p",
  "add_17p",
  "i_17q"
)

.eln2022_required_columns <- function() {
  c(
    "distinct_aberrations",
    "monosomal_karyotype",
    .eln2022_favorable,
    "t_9_11_p21_q23",
    .eln2022_adverse,
    paste0("tris", .risk_chroms),
    paste0("mono", .risk_chroms),
    names(.general_flag_patterns)
  )
}

.assign_eln2022 <- function(parsed) {
  n <- as.integer(parsed$distinct_aberrations)

  favorable <- .risk_any(parsed, .eln2022_favorable)
  kmt2a_9_11 <- .risk_any(parsed, "t_9_11_p21_q23")

  # ELN's complex definition, not the scheme-neutral complex_karyotype: the
  # hyperdiploid carve-out excludes pure-gain karyotypes. A marker counts as
  # a structural abnormality here (unlike in monosomal_karyotype, whose own
  # definition excludes it).
  n_trisomies <- .risk_count(parsed, paste0("tris", .risk_chroms))
  n_monosomies <- .risk_count(parsed, paste0("mono", .risk_chroms))
  n_structural <- .risk_count(parsed, names(.general_flag_patterns))
  hyperdiploid_gains_only <- n_trisomies >= 3L &
    n_monosomies == 0L &
    n_structural == 0L

  eln_complex <- !is.na(n) & n >= 3L & !hyperdiploid_gains_only

  adverse <- .risk_any(parsed, .eln2022_adverse) |
    eln_complex |
    .risk_any(parsed, "monosomal_karyotype")

  category <- dplyr::case_when(
    is.na(n) ~ NA_character_,
    # Class-defining, so they outrank every adverse criterion.
    favorable ~ "Favorable",
    # Must precede `adverse`: a t(9;11) also satisfies t(v;11q23).
    kmt2a_9_11 ~ "Intermediate",
    adverse ~ "Adverse",
    TRUE ~ "Intermediate"
  )

  parsed$eln2022_cyto_risk <- factor(
    category,
    levels = .eln2022_levels,
    ordered = TRUE
  )
  parsed
}
