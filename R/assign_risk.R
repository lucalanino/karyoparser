#' Assign Cytogenetic Risk Categories
#'
#' Adds a cytogenetic risk category to each row of a [parse_karyo()] result.
#' Interpretation only -- no parsing or cleaning happens here; the function
#' reads columns [parse_karyo()] already produced.
#'
#' **This is the cytogenetic component only, not a complete risk score.**
#' IPSS-R stratifies on five components -- cytogenetics, marrow blast
#' percentage, haemoglobin, platelet count, and absolute neutrophil count --
#' of which this function computes the first. A row's `ipssr_cyto_risk` is
#' therefore not an IPSS-R risk group, and must be combined with the four
#' haematological components before it means anything clinically.
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
#' @param parsed A tibble returned by [parse_karyo()]. Must carry the default
#'   [myeloid_rules] flags; if `columns` was used to narrow the output, the
#'   classes the scheme needs must still be present.
#' @param scheme Character string naming the risk scheme. Currently only
#'   `"ipssr"` (IPSS-R cytogenetic risk group).
#'
#' @return `parsed` with two columns appended:
#'   - `ipssr_cyto_risk`: ordered factor, `Very Good` < `Good` <
#'     `Intermediate` < `Poor` < `Very Poor`; `NA` for rows that cannot be
#'     scored
#'   - `ipssr_cyto_score`: integer 0-4 matching the category, `NA` alongside
#'     an `NA` category
#'
#' @examples
#' parsed <- parse_karyo(c("46,XX", "46,XX,del(5)(q13q33)", "45,XY,-7"))
#' risk <- assign_risk(parsed, scheme = "ipssr")
#' risk[, c("original_karyotype", "ipssr_cyto_risk", "ipssr_cyto_score")]
#'
#' @seealso [parse_karyo()]
#' @export
assign_risk <- function(parsed, scheme = c("ipssr")) {
  scheme <- match.arg(scheme)

  if (!is.data.frame(parsed)) {
    stop(
      "`parsed` must be a data frame returned by parse_karyo().",
      call. = FALSE
    )
  }

  required <- switch(scheme, ipssr = .ipssr_required_columns)
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

  switch(scheme, ipssr = .assign_ipssr(parsed))
}

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
  # Flags are 0/1 on parsed rows and NA on unscoreable ones; the NA rows are
  # caught by `is.na(n)` below, so flags can be flattened to 0 safely here.
  flag <- function(col) tidyr::replace_na(as.integer(parsed[[col]]), 0L)

  n <- as.integer(parsed$distinct_aberrations)
  isolated <- !is.na(n) & n == 1L
  double <- !is.na(n) & n == 2L

  abn_3q <- pmax(
    flag("inv_3_q21q26"),
    flag("inv_3_other"),
    flag("t_3q26_v"),
    flag("del_3q")
  )
  mono_7 <- flag("mono7")
  del_7q <- flag("del_7q")
  del_5q <- flag("del_5q")

  category <- dplyr::case_when(
    is.na(n) ~ NA_character_,
    n == 0L & flag("normal_karyotype") == 1L ~ "Good",
    # No scoreable token but not normal (e.g. "45,X") -- see @section above.
    n == 0L ~ NA_character_,
    n > 3L ~ "Very Poor",
    n == 3L ~ "Poor",
    double & (mono_7 == 1L | del_7q == 1L) ~ "Poor",
    double & del_5q == 1L ~ "Good",
    double ~ "Intermediate",
    isolated & mono_7 == 1L ~ "Poor",
    isolated & abn_3q == 1L ~ "Poor",
    isolated & flag("monoY") == 1L ~ "Very Good",
    isolated & flag("del_11q") == 1L ~ "Very Good",
    isolated & del_5q == 1L ~ "Good",
    isolated & flag("del_12p") == 1L ~ "Good",
    isolated & flag("del_20q") == 1L ~ "Good",
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
