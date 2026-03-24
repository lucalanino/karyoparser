# Single source of truth for all dirty-marker patterns.
# Both preprocess_karyo() and flag_unpreprocessed() loop over this list.
# Each entry:
#   detect      - character vector; fires if ANY element matches (via str_detect)
#   use_trimmed - if TRUE, detect against trimws(s) rather than raw s
#   fix         - list of list(pattern, replacement) applied in order via str_replace_all
#   detail      - human-readable description used in issue reports
.dirty_patterns <- list(
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
  leading_dotslash = list(
    detect = "^[.]+//",
    use_trimmed = TRUE,
    fix = list(
      list(pattern = "^[.]+//", replacement = "")
    ),
    detail = "String starts with './/'"
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
    detect = "\\][ .]+(?:nuc )?ish\\b",
    use_trimmed = FALSE,
    fix = list(
      list(
        pattern = "\\s+\\.?(?:nuc )?ish\\b.*$",
        replacement = ""
      )
    ),
    detail = "FISH / nuc ish suffix after last clone bracket (e.g. '[12] .nuc ish(PDGFRA x3)[20/200]')"
  ),
  trailing_narrative = list(
    detect = c("\\]\\s*\\.", "\\s+\\.\\s*[A-Z]"),
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
      # Rule 3: fallback for strings with no metaphase bracket (e.g. "46,XX .note").
      list(
        pattern = "\\s+\\.\\s*[A-Z].*$",
        replacement = ""
      )
    ),
    detail = "Trailing narrative text after last metaphase-count bracket"
  ),
  midstring_linewrap = list(
    detect = ",\\s+\\.[a-z(+]",
    use_trimmed = FALSE,
    fix = list(
      list(
        pattern = ",\\s+\\.(?=[a-z(+])",
        replacement = ","
      )
    ),
    detail = "Mid-string line-wrap artifact (e.g. ', .der(...)' or ', .+8')"
  ),
  missing_sex_comma = list(
    detect = ",(XXXXY|XXXX|XXX|XXY|XYY|XY|XX|X|Y)\\s+(?=[a-z(+])",
    use_trimmed = FALSE,
    fix = list(
      list(
        pattern = "(,(XXXXY|XXXX|XXX|XXY|XYY|XY|XX|X|Y))\\s+(?=[a-z(+])",
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
  )
)

#' Clean Dirty ISCN Karyotype Strings
#'
#' Strips clinical free text, database artifacts, and HTML entities from raw
#' karyotype strings before parsing. Runs BEFORE `normalize_iscn()`, which
#' handles unicode/notation normalization inside `parse_karyo()`.
#'
#' Rules applied in order:
#' 1. Trim leading/trailing whitespace and collapse newlines/tabs
#' 2. Decode HTML entities (`&lt;` → `<`, `&gt;` → `>`, `&amp;` → `&`)
#' 3. Strip leading `.//` (and variants like `..//`)
#' 4. Strip leading dot(s) before a digit (e.g. `.46,XX` → `46,XX`)
#' 5. Strip FISH/nuc ish suffix after last clone bracket (e.g. `[12] .nuc ish(...)`)
#' 6. Strip trailing narrative: `] .text` → `]`; space-dot-capital mid-string
#' 7. Collapse mid-string line-wrap artifacts (`, .der(...)` → `,der(...)`; `, .+8` → `,+8`)
#' 8. Insert missing comma after sex chromosome complement (`46,XX der(...)` → `46,XX,der(...)`)
#' 9. Remove space between count and `mar` token (`+1~4 mar` → `+1~4mar`)
#' 10. Trim again
#'
#' @param x Character vector of raw karyotype strings.
#' @return Character vector of cleaned karyotype strings (same length as `x`).
#' @export
preprocess_karyo <- function(x) {
  if (length(x) == 0) {
    return(x)
  }
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "[\n\r\t]+", " ")
  for (nm in names(.dirty_patterns)) {
    dp <- .dirty_patterns[[nm]]
    for (fx in dp$fix) {
      x <- stringr::str_replace_all(x, fx$pattern, fx$replacement)
    }
  }
  x <- stringr::str_trim(x)
  x
}

#' Flag Unpreprocessed ISCN Strings
#'
#' Detects dirty markers in karyotype strings that `normalize_iscn()` cannot
#' fix. Called internally by `check_karyo()`.
#'
#' Detected issue types:
#' - `html_entities`: contains `&lt;`, `&gt;`, or `&amp;`
#' - `leading_dotslash`: string starts with `.//` (after trimming)
#' - `leading_dot`: string starts with dot(s) before a digit
#' - `fish_notation`: FISH/nuc ish suffix after last clone bracket
#' - `trailing_narrative`: bracket followed by dot (`] .text`) OR
#'   space-dot-capital pattern (` .Text`, no bracket required)
#' - `midstring_linewrap`: `, .lowercase` or `, .+` mid-string line-wrap artifact
#' - `missing_sex_comma`: sex complement followed by space instead of comma
#'   (e.g. `46,XX der(...)` → `46,XX,der(...)`)
#' - `mar_space`: space between count and `mar` token (e.g. `+1~4 mar`)
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
    return(tibble::tibble(
      row_index = integer(),
      karyotype = character(),
      issue_type = character(),
      issue_detail = character()
    ))
  }
  dplyr::bind_rows(issue_list[seq_len(n_issues)])
}

#' Apply preprocess_karyo() to a subset of rows
#'
#' @param vec Character vector.
#' @param row_indices Integer indices of elements to clean.
#' @return `vec` with `row_indices` elements replaced by cleaned versions.
#' @keywords internal
apply_preprocess_to_rows <- function(vec, row_indices) {
  vec[row_indices] <- preprocess_karyo(vec[row_indices])
  vec
}

#' Preprocess ISCN Karyotype Strings
#'
#' Normalizes unicode characters, whitespace, delimiters, and notation
#' variants (idem, sl, psu dic) in ISCN karyotype strings. Runs automatically
#' inside `parse_karyo()`.
#'
#' @param x Character vector of raw karyotype strings.
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
