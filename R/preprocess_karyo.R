#' Clean and Normalize ISCN Karyotype Strings
#'
#' The complete cleaning and normalization pipeline. Fixes dirty markers then
#' runs `normalize_iscn()` as a final pass. Output is fully normalized and
#' parser-ready. This is the only place in the workflow where cleaning or
#' normalization occurs -- `parse_karyo()` does none.
#'
#' Rules applied in order:
#' 1. Trim leading/trailing whitespace
#' 2. Replace unicode spaces (NBSP), dashes (em/en-dash), fullwidth characters
#' 3. Collapse embedded newlines/tabs to spaces
#' 4. Decode HTML entities (`&lt;` -> `<`, `&gt;` -> `>`, `&amp;` -> `&`)
#' 5. Strip leading dot(s) before a digit (e.g. `.46,XX` -> `46,XX`)
#' 6. Insert missing/dotted separator between chromosome count and sex
#'    complement (e.g. `46XY` or `45.XY` -> `46,XY`/`45,XY`)
#' 7. Strip FISH/nuc ish annotation (suffix or mid-clone before metaphase count)
#' 8. Strip trailing narrative: `] .text` -> `]`; `) Capital text` -> `)`
#' 9. Collapse mid-string line-wrap artifacts (`, .der(...)` -> `,der(...)`)
#' 10. Insert missing comma after sex chromosome complement
#' 11. Remove space between count and `mar` token (`+1~4 mar` -> `+1~4mar`)
#' 12. `normalize_iscn()`: whitespace collapsing, delimiter tightening,
#'     idem/sl/cp normalization
#'
#' Note: `zero_host_chimera` strings (`.//` or `//` prefix) are donor-only
#' chimeras with no host metaphases. The leading `.//` prefix is stripped,
#' leaving the donor clone(s), which are then parsed (under `on_chimeric =
#' "host"` these rows are instead returned NA).
#'
#' @param karyotypes Character vector of raw karyotype strings, the
#'   `karyo_check` tibble returned by `check_karyo()`, or a data frame
#'   containing a karyotype
#'   column. When a `karyo_check` tibble is supplied, the assessment it already
#'   computed is reused directly -- no re-scanning -- and any id column detected
#'   upstream is propagated automatically. When a plain data frame is supplied,
#'   the karyotype column is auto-detected or specified via `karyotype_column`.
#' @param karyotype_column Character. Name of the karyotype column when `x` is
#'   a plain data frame. If `NULL` (default), auto-detected from common names.
#'   Ignored for character vector or `karyo_check` input.
#' @param id_column Character. Name of the id column when `x` is a plain data
#'   frame. If `NULL` (default), auto-detected from common names and used
#'   silently. Ignored for character vector or `karyo_check` input (id is
#'   propagated automatically from `karyo_check` attrs in that case).
#' @param on_chimeric Character string controlling which clone is kept for
#'   chimeric karyotypes (those containing a `//` separator). One of:
#'   - `"default"`: host clone (before `//`) for normal chimeras; donor clone
#'     (after `//`) for zero-host chimeras (`.//` prefix).
#'   - `"host"`: host clone only; zero-host chimeras become NA (the NA is
#'     policy-induced, so the row stays classified as fixable, not unfixable).
#'   - `"donor"`: everything after `//` (the donor population).
#'   Rows with two or more `//` separators are always unfixable. The choice is
#'   baked into the `preprocessed` column and recorded so `parse_karyo()` can
#'   reuse it.
#' @param verbose Logical. If `TRUE`, prints a processing summary. Default
#'   `FALSE`.
#' @return A `karyo_preprocessed` tibble. Columns: optional id column (first,
#'   when detected or propagated), `original` (raw input, always unchanged),
#'   `preprocessed` (fully normalized string, `NA_character_` for unfixable
#'   rows), `status` (`"clean"`, `"fixed"`, or `"unfixable"`).
#'
#'   The returned tibble can be passed directly to `parse_karyo()` without
#'   specifying `karyotype_column` or `id_column` -- both are inferred
#'   automatically from the object's class and cached attributes.
#' @export
preprocess_karyo <- function(
  karyotypes,
  karyotype_column = NULL,
  id_column = NULL,
  on_chimeric = c("default", "host", "donor"),
  verbose = FALSE
) {
  on_chimeric <- match.arg(on_chimeric)
  id_col_name <- NULL
  id_values <- NULL
  if (inherits(karyotypes, "karyo_check")) {
    cached_assessment <- attr(karyotypes, ".kp_assessment")
    id_col_name <- attr(karyotypes, ".kp_id_col")
    id_values <- attr(karyotypes, ".kp_id_values")
    karyotypes <- attr(karyotypes, ".kp_raw")
  } else if (is.data.frame(karyotypes)) {
    extracted <- .extract_df_input(
      karyotypes,
      karyotype_column,
      id_column,
      verbose,
      "preprocess_karyo"
    )
    id_col_name <- extracted$id_col_name
    id_values <- extracted$id_values
    karyotypes <- extracted$raw_vec
    cached_assessment <- NULL
  } else {
    if (!is.character(karyotypes)) {
      stop(
        "`karyotypes` must be a character vector, data frame, or karyo_check object.",
        call. = FALSE
      )
    }
    cached_assessment <- NULL
  }

  if (length(karyotypes) == 0) {
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
    attr(out, ".kp_chimeric_clone") <- character(0)
    attr(out, ".kp_on_chimeric") <- on_chimeric
    attr(out, ".kp_id_col") <- id_col_name
    attr(out, ".kp_id_values") <- id_values
    return(out)
  }

  n <- length(karyotypes)
  if (isTRUE(verbose)) {
    message("Preprocessing ", n, " karyotype(s)...")
  }

  # The cached assessment from check_karyo() always uses the default clone
  # selection; recompute if a different on_chimeric policy is requested.
  assessment <- if (
    !is.null(cached_assessment) &&
      identical(cached_assessment$on_chimeric, on_chimeric)
  ) {
    cached_assessment
  } else {
    .assess_karyotypes(karyotypes, on_chimeric)
  }
  processed <- assessment$processed

  # Status: unfixable > fixed (dirty or chimeric, successfully resolved) > clean
  fixed_row_indices <- setdiff(
    unique(c(
      assessment$dirty_row_indices,
      assessment$chimeric_row_indices,
      assessment$zhc_row_indices
    )),
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
    original = as.character(karyotypes),
    preprocessed = as.character(processed),
    status = status
  )

  if (!is.null(id_col_name)) {
    out[[id_col_name]] <- id_values
    out <- dplyr::relocate(out, dplyr::all_of(id_col_name), .before = 1)
  }

  # Tag output so parse_karyo() can skip .assess_karyotypes() on the result.
  class(out) <- c("karyo_preprocessed", class(out))
  attr(out, ".kp_fixable_rows") <- unique(c(
    assessment$dirty_row_indices,
    assessment$chimeric_row_indices,
    assessment$zhc_row_indices
  ))
  attr(out, ".kp_unfixable_rows") <- assessment$unfixable_row_indices
  attr(out, ".kp_chimeric_rows") <- assessment$chimeric_row_indices
  attr(out, ".kp_chimeric_all_rows") <- unique(c(
    assessment$chimeric_row_indices,
    assessment$zhc_row_indices,
    assessment$multi_row_indices
  ))
  attr(out, ".kp_chimeric_clone") <- assessment$chimeric_clone
  attr(out, ".kp_on_chimeric") <- on_chimeric
  attr(out, ".kp_id_col") <- id_col_name
  attr(out, ".kp_id_values") <- id_values
  out
}

# Final normalization: whitespace, delimiters, idem/sl/cp/range variants.
normalize_iscn <- function(x) {
  if (length(x) == 0) {
    return(x)
  }
  x <- stringr::str_replace_all(x, "[\u00A0\u2007\u202F]", " ")
  x <- stringr::str_replace_all(
    x,
    "[\u2212\u2012\u2013\u2014\uFE63\uFF0D]",
    "-"
  ) # minus/dashes -> '-'
  x <- stringr::str_replace_all(x, "[\uFF0B]", "+")
  x <- stringr::str_replace_all(x, "\n|\r|\t", " ")
  x <- stringr::str_replace_all(x, " +", " ")
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "\\s*,\\s*", ",")
  x <- stringr::str_replace_all(x, "\\s*/\\s*", "/")
  x <- stringr::str_replace_all(x, "\\s*;\\s*", ";")
  x <- stringr::str_replace_all(x, "\\s*\\)\\s*", ")")
  x <- stringr::str_replace_all(x, "\\s*\\(\\s*", "(")
  x <- stringr::str_replace_all(x, "\\s+\\[", "[")
  x <- stringr::str_replace_all(x, "(?i)\\bpsu\\s*dic\\b", "psu dic")
  x <- stringr::str_replace_all(x, "\\s*\\bcp\\s*\\[(\\d+)\\]", "[cp\\1]")
  x <- stringr::str_replace_all(x, "\\[cp\\s*(\\d+)\\]", "[cp\\1]")
  x <- stringr::str_replace_all(x, "^(\\d+)-(\\d+)(?=,)", "\\1~\\2")
  x <- stringr::str_replace_all(x, "(/)(\\d+)-(\\d+)(?=,)", "\\1\\2~\\3")
  x <- stringr::str_replace_all(x, "(?i)\\bIDEM\\b", "idem")
  x <- stringr::str_replace_all(x, "(?i)\\bSL\\b", "sl")
  x <- stringr::str_replace_all(x, ",{2,}", ",")
  x <- stringr::str_replace_all(x, "^[,;/\\.]+", "")
  x <- stringr::str_replace_all(x, "[,;/\\.]+$", "")
  x
}
