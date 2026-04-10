# Sex complement alternation derived from the canonical constant, longest variants
# first so the regex engine doesn't match a shorter prefix before a longer one.
.sex_alt <- paste(
  .sex_complements[order(-nchar(.sex_complements))],
  collapse = "|"
)

# Single source of truth for all pre-normalization issue patterns.
# Both preprocess_karyo() and flag_unpreprocessed() loop over this list.
# Each entry:
#   detect      - character vector; fires if ANY element matches (via str_detect)
#   use_trimmed - if TRUE, detect against trimws(s) rather than raw s
#   fix         - list of list(pattern, replacement) applied in order via str_replace_all
#               - empty list means the issue is detected but NOT fixable by preprocess_karyo()
#   detail      - human-readable description used in issue reports
.dirty_patterns <- list(
  unicode_notation = list(
    detect = "[\u00A0\u2007\u202F\u2212\u2012\u2013\u2014\uFE63\uFF0D\uFF0B]",
    use_trimmed = FALSE,
    fix = list(
      list(pattern = "[\u00A0\u2007\u202F]", replacement = " "),
      list(
        pattern = "[\u2212\u2012\u2013\u2014\uFE63\uFF0D]",
        replacement = "-"
      ),
      list(pattern = "\uFF0B", replacement = "+")
    ),
    detail = "Contains unicode spaces (NBSP), dashes (em/en-dash), or fullwidth characters"
  ),
  embedded_newline = list(
    detect = "[\n\r\t]",
    use_trimmed = FALSE,
    fix = list(
      list(pattern = "[\n\r\t]+", replacement = " ")
    ),
    detail = "Contains embedded newline or tab characters"
  ),
  html_entities = list(
    detect = "&lt;|&gt;|&amp;",
    use_trimmed = FALSE,
    fix = list(
      list(pattern = "&lt;", replacement = "<"),
      list(pattern = "&gt;", replacement = ">"),
      list(pattern = "&amp;", replacement = "&")
    ),
    detail = "Contains HTML entities (&lt;, &gt;, or &amp;)"
  ),
  leading_dot = list(
    detect = "^[.]+(?=\\d)",
    use_trimmed = TRUE,
    fix = list(
      list(pattern = "^[.]+(?=\\d)", replacement = "")
    ),
    detail = "String starts with dot(s) before chromosome count"
  ),
  fish_notation = list(
    detect = "\\][. ]+(?:nuc )?ish\\b",
    use_trimmed = FALSE,
    fix = list(
      list(
        pattern = "[. ]+(?:nuc )?ish\\b.*$",
        replacement = ""
      )
    ),
    detail = "FISH / nuc ish suffix after last clone bracket (e.g. '[12] .nuc ish(PDGFRA x3)[20/200]' or '[1].ish t(...)')"
  ),
  trailing_narrative = list(
    detect = c("\\]\\s*\\.", "\\s+\\.\\s*[A-Z]", "\\)\\s+[A-Z][a-z]"),
    use_trimmed = FALSE,
    fix = list(
      # Rule 1: strip everything after the last metaphase-count bracket.
      # Greedy .* anchors to the rightmost [n], [cpN], or [n~m] bracket.
      list(
        pattern = "^(.*\\[(?:cp)?\\d+(?:[~-]\\d+)?\\]).*$",
        replacement = "\\1"
      ),
      # Rule 2: collapse ] ./ separator artifact left between clones.
      list(
        pattern = "\\]\\s*\\./",
        replacement = "]/"
      ),
      # Rule 3: strip narrative after last ')' when no metaphase bracket present
      # (e.g. "46,XX,t(9;22)(q34;q11.2) Abnormal female karyotype").
      # Requires Capital+lowercase to avoid false positives on ISCN tokens.
      list(
        pattern = "^(.*\\))\\s+[A-Z][a-z].*$",
        replacement = "\\1"
      ),
      # Rule 4: fallback for strings with no metaphase bracket or parens.
      list(
        pattern = "\\s+\\.\\s*[A-Z].*$",
        replacement = ""
      )
    ),
    detail = "Trailing narrative text after last metaphase-count bracket or closing parenthesis"
  ),
  midstring_linewrap = list(
    detect = ",\\s+\\.[a-z(+]",
    use_trimmed = FALSE,
    fix = list(
      list(
        pattern = ",\\s+\\.?(?=[a-z(+])",
        replacement = ","
      )
    ),
    detail = "Mid-string line-wrap artifact (e.g. ', .der(...)' or bare ', +8' without dot)"
  ),
  missing_sex_comma = list(
    detect = paste0(",(", .sex_alt, ")\\s+(?=[a-z(+])"),
    use_trimmed = FALSE,
    fix = list(
      list(
        pattern = paste0("(,(", .sex_alt, "))\\s+(?=[a-z(+])"),
        replacement = "\\1,"
      )
    ),
    detail = "Missing comma between sex chromosome complement and first aberration (e.g. '46,XX der(...)' should be '46,XX,der(...)')"
  ),
  mar_space = list(
    detect = "[+~0-9-] mar\\b",
    use_trimmed = FALSE,
    fix = list(
      list(
        pattern = "([+~0-9-]) (mar\\b)",
        replacement = "\\1\\2"
      )
    ),
    detail = "Space between count and 'mar' token (e.g. '+1~4 mar' should be '+1~4mar')"
  ),
  # Detected in flag_unpreprocessed() on the RAW string, before any fixes or
  # normalize_iscn() runs. normalize_iscn() strips leading punctuation
  # (^[,;/.]+) which would destroy the signal.
  zero_host_chimera = list(
    detect = "^[.]*//",
    use_trimmed = TRUE,
    fix = list(),
    detail = "String starts with './/'' or '//': donor-only chimera with no host metaphases"
  )
)

# Issue types resolvable by parse_karyo() under on_issues="fix".
# Defined here (after .dirty_patterns) so load order is guaranteed.
.fixable_issue_types <- c(
  names(Filter(function(p) length(p$fix) > 0, .dirty_patterns)),
  "chimeric_separator"
)

# Complete fixed schema of all detectable issue types, in display order.
.all_issue_types <- c(
  names(.dirty_patterns),
  "empty",
  "no_chromosome_count",
  "chimeric_separator",
  "updated_iscn",
  "unbalanced_parentheses",
  "unbalanced_brackets",
  "no_sex_complement",
  "invalid_idem",
  "unparseable_bracket"
)

#' Clean and Normalize ISCN Karyotype Strings
#'
#' The complete cleaning and normalization pipeline. Fixes dirty markers then
#' runs `normalize_iscn()` as a final pass. Output is fully normalized and
#' parser-ready. This is the only place in the workflow where cleaning or
#' normalization occurs — `parse_karyo()` does none.
#'
#' Rules applied in order:
#' 1. Trim leading/trailing whitespace
#' 2. Replace unicode spaces (NBSP), dashes (em/en-dash), fullwidth characters
#' 3. Collapse embedded newlines/tabs to spaces
#' 4. Decode HTML entities (`&lt;` → `<`, `&gt;` → `>`, `&amp;` → `&`)
#' 5. Strip leading dot(s) before a digit (e.g. `.46,XX` → `46,XX`)
#' 6. Strip FISH/nuc ish suffix after last clone bracket
#' 7. Strip trailing narrative: `] .text` → `]`; `) Capital text` → `)`
#' 8. Collapse mid-string line-wrap artifacts (`, .der(...)` → `,der(...)`)
#' 9. Insert missing comma after sex chromosome complement
#' 10. Remove space between count and `mar` token (`+1~4 mar` → `+1~4mar`)
#' 11. `normalize_iscn()`: whitespace collapsing, delimiter tightening,
#'     idem/sl/cp normalization
#'
#' Note: `zero_host_chimera` strings (`.//` or `//` prefix) are detected but
#' not modified — their `fix` list is empty, so the loop skips them. Always NA.
#'
#' @param x Character vector of raw karyotype strings, the `karyo_check`
#'   tibble returned by `check_karyo()`, or a data frame containing a karyotype
#'   column. When a `karyo_check` tibble is supplied, the assessment it already
#'   computed is reused directly — no re-scanning — and any id column detected
#'   upstream is propagated automatically. When a plain data frame is supplied,
#'   the karyotype column is auto-detected or specified via `karyotype_column`.
#' @param karyotype_column Character. Name of the karyotype column when `x` is
#'   a plain data frame. Ignored for character vector or `karyo_check` input.
#' @param id_column Character. Name of the id column when `x` is a plain data
#'   frame. Ignored for character vector or `karyo_check` input (id is
#'   propagated automatically from `karyo_check` attrs in that case).
#' @param verbose Logical. If `TRUE`, prints a processing summary. Default
#'   `FALSE`.
#' @return A `karyo_preprocessed` tibble. Columns: optional id column (first,
#'   when detected or propagated), `original` (raw input, always unchanged),
#'   `preprocessed` (fully normalized string, `NA_character_` for unfixable
#'   rows), `status` (`"clean"`, `"fixed"`, or `"unfixable"`).
#'
#'   The returned tibble can be passed directly to `parse_karyo()` without
#'   specifying `karyotype_column` or `id_column` — both are inferred
#'   automatically from the object's class and cached attributes.
#' @export
preprocess_karyo <- function(
  x,
  karyotype_column = NULL,
  id_column = NULL,
  verbose = FALSE
) {
  # Input routing
  id_col_name <- NULL
  id_values <- NULL
  if (inherits(x, "karyo_check")) {
    # Reuse cached assessment and propagate id from upstream check_karyo()
    cached_assessment <- attr(x, ".kp_assessment")
    id_col_name <- attr(x, ".kp_id_col")
    id_values <- attr(x, ".kp_id_values")
    x <- attr(x, ".kp_raw")
  } else if (is.data.frame(x)) {
    extracted <- .extract_df_input(
      x,
      karyotype_column,
      id_column,
      verbose,
      "preprocess_karyo"
    )
    id_col_name <- extracted$id_col_name
    id_values <- extracted$id_values
    x <- extracted$raw_vec
    cached_assessment <- NULL
  } else {
    if (!is.character(x)) {
      stop(
        "`x` must be a character vector, data frame, or karyo_check object.",
        call. = FALSE
      )
    }
    cached_assessment <- NULL
  }

  if (length(x) == 0) {
    out <- tibble::tibble(
      original = character(),
      preprocessed = character(),
      status = character()
    )
    if (!is.null(id_col_name)) {
      out[[id_col_name]] <- character()
      out <- dplyr::relocate(out, dplyr::all_of(id_col_name), .before = 1)
    }
    class(out) <- c("karyo_preprocessed", class(out))
    attr(out, ".kp_fixable_rows") <- integer(0)
    attr(out, ".kp_unfixable_rows") <- integer(0)
    attr(out, ".kp_chimeric_rows") <- integer(0)
    attr(out, ".kp_chimeric_all_rows") <- integer(0)
    attr(out, ".kp_id_col") <- id_col_name
    attr(out, ".kp_id_values") <- id_values
    return(out)
  }

  n <- length(x)
  if (isTRUE(verbose)) {
    message("Preprocessing ", n, " karyotype(s)...")
  }

  assessment <- if (!is.null(cached_assessment)) {
    cached_assessment
  } else {
    .assess_karyotypes(x)
  }
  processed <- assessment$processed

  # Status: unfixable > fixed (dirty or chimeric, successfully resolved) > clean
  fixed_row_indices <- setdiff(
    union(assessment$dirty_row_indices, assessment$chimeric_row_indices),
    assessment$unfixable_row_indices
  )
  status <- ifelse(
    seq_len(n) %in% assessment$unfixable_row_indices,
    "unfixable",
    ifelse(seq_len(n) %in% fixed_row_indices, "fixed", "clean")
  )

  if (isTRUE(verbose)) {
    n_fixed <- sum(status == "fixed")
    n_unfixable <- length(assessment$unfixable_row_indices)
    if (n_fixed == 0 && n_unfixable == 0) {
      message("  All clean.")
    } else {
      n_clean <- sum(status == "clean", na.rm = TRUE)
      if (n_fixed > 0) {
        message("  Fixed:      ", n_fixed)
      }
      if (n_clean > 0) {
        message("  Clean:      ", n_clean)
      }
      if (n_unfixable > 0) {
        message(
          "  Unfixable:  ",
          n_unfixable,
          "  \u2014 run check_karyo() to investigate"
        )
      }
    }
  }

  out <- tibble::tibble(
    original = as.character(x),
    preprocessed = as.character(processed),
    status = status
  )

  # Add id column as first column when present
  if (!is.null(id_col_name)) {
    out[[id_col_name]] <- id_values
    out <- dplyr::relocate(out, dplyr::all_of(id_col_name), .before = 1)
  }

  # Tag output so parse_karyo() can skip .assess_karyotypes() on the result.
  zhc_row_indices <- unique(
    assessment$reported_issues$row_index[
      assessment$reported_issues$issue_type == "zero_host_chimera"
    ]
  )
  class(out) <- c("karyo_preprocessed", class(out))
  attr(out, ".kp_fixable_rows") <- unique(
    union(assessment$dirty_row_indices, assessment$chimeric_row_indices)
  )
  attr(out, ".kp_unfixable_rows") <- assessment$unfixable_row_indices
  attr(out, ".kp_chimeric_rows") <- assessment$chimeric_row_indices
  attr(out, ".kp_chimeric_all_rows") <- unique(
    c(assessment$chimeric_row_indices, zhc_row_indices)
  )
  attr(out, ".kp_id_col") <- id_col_name
  attr(out, ".kp_id_values") <- id_values
  out
}

#' Flag Unpreprocessed ISCN Strings
#'
#' Detects pre-normalization issues in raw karyotype strings: dirty markers
#' that `preprocess_karyo()` can fix, and structural issues that must be caught
#' before `normalize_iscn()` destroys the signal. Called internally by
#' `check_karyo()`.
#'
#' Detected issue types:
#' - `unicode_notation`: unicode spaces (NBSP), dashes (em/en-dash), or
#'   fullwidth characters
#' - `embedded_newline`: embedded `\n`, `\r`, or `\t` characters
#' - `html_entities`: contains `&lt;`, `&gt;`, or `&amp;`
#' - `leading_dot`: string starts with dot(s) before a digit
#' - `fish_notation`: FISH/nuc ish suffix after last clone bracket
#' - `trailing_narrative`: bracket followed by dot (`] .text`) OR
#'   space-dot-capital pattern (` .Text`, no bracket required)
#' - `midstring_linewrap`: `, .lowercase` or `, .+` mid-string line-wrap artifact
#' - `missing_sex_comma`: sex complement followed by space instead of comma
#'   (e.g. `46,XX der(...)` → `46,XX,der(...)`)
#' - `mar_space`: space between count and `mar` token (e.g. `+1~4 mar`)
#' - `zero_host_chimera`: string starts with `.//` or `//` — donor-only chimera
#'   with no host metaphases; has `fix = list()` so always unfixable/NA
#'
#' @param x Character vector of karyotype strings.
#' @return A tibble with columns: `row_index`, `karyotype` (truncated to 40
#'   chars), `issue_type`, `issue_detail`. Each issue generates one row;
#'   multiple issues on the same input row appear as separate rows.
#' @keywords internal
flag_unpreprocessed <- function(x) {
  issue_list <- vector("list", length(x) * length(.dirty_patterns))
  n_issues <- 0L

  for (i in seq_along(x)) {
    s <- x[i]
    if (is.na(s)) {
      next
    }
    s_trim <- trimws(s)

    for (nm in names(.dirty_patterns)) {
      dp <- .dirty_patterns[[nm]]
      target <- if (isTRUE(dp$use_trimmed)) s_trim else s
      if (any(stringr::str_detect(target, dp$detect))) {
        n_issues <- n_issues + 1L
        issue_list[[n_issues]] <- tibble::tibble(
          row_index = i,
          karyotype = truncate_str(s),
          issue_type = nm,
          issue_detail = dp$detail
        )
      }
    }
  }

  if (n_issues == 0L) {
    return(empty_issues_tibble())
  }
  dplyr::bind_rows(issue_list[seq_len(n_issues)])
}

#' Normalize ISCN Notation
#'
#' Final normalization pass: collapses whitespace, tightens delimiters,
#' normalizes notation variants (idem, sl, cp, psu dic, range notation).
#' Called at the end of `preprocess_karyo()` and internally by
#' `.assess_karyotypes()` for structural validation checks.
#' Never called directly by `parse_karyo()`.
#'
#' @param x Character vector of karyotype strings.
#' @return Character vector of normalized karyotype strings.
#' @keywords internal
normalize_iscn <- function(x) {
  if (length(x) == 0) {
    return(x)
  }
  # normalize unicode spaces and signs
  x <- stringr::str_replace_all(x, "[\u00A0\u2007\u202F]", " ") # NBSPs
  x <- stringr::str_replace_all(
    x,
    "[\u2212\u2012\u2013\u2014\uFE63\uFF0D]",
    "-"
  ) # minus/dashes -> '-'
  x <- stringr::str_replace_all(x, "[\uFF0B]", "+") # fullwidth '+' -> '+'
  # minimal whitespace policy
  x <- stringr::str_replace_all(x, "\n|\r|\t", " ")
  x <- stringr::str_replace_all(x, " +", " ")
  x <- stringr::str_trim(x)
  # tighten around delimiters
  x <- stringr::str_replace_all(x, "\\s*,\\s*", ",")
  x <- stringr::str_replace_all(x, "\\s*/\\s*", "/")
  x <- stringr::str_replace_all(x, "\\s*;\\s*", ";")
  x <- stringr::str_replace_all(x, "\\s*\\)\\s*", ")")
  x <- stringr::str_replace_all(x, "\\s*\\(\\s*", "(")
  x <- stringr::str_replace_all(x, "\\s+\\[", "[")
  # normalizations
  x <- stringr::str_replace_all(x, "(?i)\\bpsu\\s*dic\\b", "psu dic")
  x <- stringr::str_replace_all(x, "\\s*\\bcp\\s*\\[(\\d+)\\]", "[cp\\1]")
  x <- stringr::str_replace_all(x, "\\[cp\\s*(\\d+)\\]", "[cp\\1]")
  x <- stringr::str_replace_all(x, "^(\\d+)-(\\d+)(?=,)", "\\1~\\2")
  x <- stringr::str_replace_all(x, "(/)(\\d+)-(\\d+)(?=,)", "\\1\\2~\\3")
  x <- stringr::str_replace_all(x, "(?i)\\bIDEM\\b", "idem")
  x <- stringr::str_replace_all(x, "(?i)\\bSL\\b", "sl")
  # collapse duplicate commas; drop leading/trailing delimiters
  x <- stringr::str_replace_all(x, ",{2,}", ",")
  x <- stringr::str_replace_all(x, "^[,;/\\.]+", "")
  x <- stringr::str_replace_all(x, "[,;/\\.]+$", "")
  x
}
