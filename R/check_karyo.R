#' Check Karyotype Strings for Issues
#'
#' Checks raw karyotype strings for formatting artifacts and structural errors.
#' Read-only: this is a diagnostic step and never modifies a karyotype string
#' (cleaning happens in `preprocess_karyo()`). Returns a wide-format tibble
#' with one row per input string: each possible issue type is a column
#' (`0`/`1`), plus summary `fixable` and `unfixable` columns. When `verbose =
#' TRUE`, prints a count summary and per-issue-type breakdown.
#'
#' Fixable issues (resolvable by `parse_karyo()` under
#' `on_issues = "preprocess"`):
#' dirty markers (`unicode_notation`, `embedded_newline`, `html_entities`,
#' `leading_dot`, `count_sex_separator`, `fish_notation`, `trailing_narrative`,
#' `midstring_linewrap`, `missing_sex_comma`, `mar_space`, `non_clonal_sca`),
#' `chimeric_separator`, and `zero_host_chimera` (the donor clone is parsed;
#' see `on_chimeric` in `parse_karyo()`).
#'
#' Unfixable issues (always returned as NA by `parse_karyo()`):
#' `multiple_chimeric_separator`, `empty`, `no_chromosome_count`,
#' `updated_iscn`, `unbalanced_parentheses`, `unbalanced_brackets`,
#' `no_sex_complement`, `single_token`, `constitutional_sex_complement`,
#' `mosaic_karyotype`, `invalid_idem`, `unparseable_bracket`.
#'
#' `single_token` catches a bare comma-less string (e.g. `"8"`): it's
#' ambiguous whether that was meant as a chromosome count or as an
#' abnormality whose leading `+`/`-` sign was stripped (a common artifact of
#' opening ISCN strings in Excel, which silently drops a redundant leading
#' `+` from numeric-looking cells).
#'
#' Some constructs are deliberately out of scope and flagged unfixable rather
#' than parsed: constitutional abnormalities (a sex complement with a `c`
#' suffix, e.g. `47,XXYc`, flagged `constitutional_sex_complement`) and mosaic
#' karyotypes (a leading `mos` prefix, flagged `mosaic_karyotype`). These take
#' precedence over `no_chromosome_count` and `no_sex_complement` so a present
#' count or complement is not mislabelled.
#'
#' Non-clonal single-cell abnormalities (an `ncSCA` token, standalone or
#' parenthesized, flagged `non_clonal_sca`) carry no structural detail of
#' their own, so the token is silently stripped and the remaining clone(s)
#' parsed as usual.
#'
#' @param karyotypes Character vector of karyotype strings, or a data frame
#'   containing a karyotype column. If a data frame, the karyotype column is
#'   auto-detected from common names (`karyotype`, `iscn`, etc.) or
#'   specified via `karyotype_column`. An id column is also auto-detected
#'   or specified via
#'   `id_column`; when found it is included as the first column of the output
#'   and propagated through subsequent pipeline steps.
#' @param karyotype_column Character. Name of the karyotype column when
#'   `karyotypes` is a data frame. If `NULL` (default), auto-detected from
#'   common names. Ignored when `karyotypes` is a character vector.
#' @param id_column Character. Name of the id column when `karyotypes` is a
#'   data frame. If `NULL` (default), auto-detected from common names
#'   (`sample_id`, `patient_id`, `id`, `mrn`, etc.) and used silently. Ignored
#'   when `karyotypes` is a character vector.
#' @param verbose Logical. If `TRUE`, prints a count summary and per-issue-type
#'   breakdown. Default `FALSE`.
#' @return A `karyo_check` tibble with `length(karyotypes)` rows (or
#'   `nrow(karyotypes)` when input is a data frame). Columns: optional id column
#'   (first, when detected), `karyotype` (full input string), `fixable`,
#'   `unfixable`, then one integer column per unfixable issue type
#'   (alphabetical), then one integer column per fixable issue type
#'   (alphabetical). The tibble can be passed directly to
#'   `preprocess_karyo()`, which will reuse the cached assessment and
#'   propagate the id column without re-scanning.
#' @export
check_karyo <- function(
  karyotypes,
  karyotype_column = NULL,
  id_column = NULL,
  verbose = FALSE
) {
  if (inherits(karyotypes, "karyo_check")) {
    stop("`karyotypes` is already a karyo_check object.", call. = FALSE)
  }

  id_col_name <- NULL
  id_values <- NULL
  karyotype_col_name <- NULL
  if (is.data.frame(karyotypes)) {
    extracted <- .extract_df_input(
      karyotypes,
      karyotype_column,
      id_column,
      verbose,
      "check_karyo"
    )
    id_col_name <- extracted$id_col_name
    id_values <- extracted$id_values
    karyotype_col_name <- extracted$karyotype_col_name
    karyotypes <- extracted$raw_vec
  } else if (!is.character(karyotypes)) {
    stop(
      "`karyotypes` must be a character vector or data frame.",
      call. = FALSE
    )
  }

  n <- length(karyotypes)
  unfixable_issue_cols <- sort(setdiff(.all_issue_types, .fixable_issue_types))
  fixable_issue_cols <- sort(intersect(.fixable_issue_types, .all_issue_types))
  all_cols <- c(
    "karyotype",
    "fixable",
    "unfixable",
    unfixable_issue_cols,
    fixable_issue_cols
  )

  if (n == 0) {
    out <- tibble::tibble(karyotype = character())
    for (nm in c(
      "fixable",
      "unfixable",
      unfixable_issue_cols,
      fixable_issue_cols
    )) {
      out[[nm]] <- integer()
    }
    out <- out[, all_cols]
    if (!is.null(id_col_name)) {
      out[[id_col_name]] <- character()
      out <- dplyr::relocate(out, dplyr::all_of(id_col_name), .before = 1)
    }
    class(out) <- c("karyo_check", class(out))
    attr(out, ".kp_assessment") <- NULL
    attr(out, ".kp_raw") <- karyotypes
    attr(out, ".kp_karyotype_col") <- karyotype_col_name
    attr(out, ".kp_id_col") <- id_col_name
    attr(out, ".kp_id_values") <- id_values
    return(out)
  }

  if (isTRUE(verbose)) {
    message("Checking ", n, " karyotype(s)...")
  }

  assessment <- .assess_karyotypes(karyotypes)
  long <- assessment$reported_issues

  out <- tibble::tibble(karyotype = as.character(karyotypes))
  for (nm in c(unfixable_issue_cols, fixable_issue_cols)) {
    rows_with_type <- long$row_index[long$issue_type == nm]
    out[[nm]] <- as.integer(seq_len(n) %in% rows_with_type)
  }

  out$fixable <- as.integer(
    rowSums(out[, fixable_issue_cols, drop = FALSE]) > 0 &
      !seq_len(n) %in% assessment$unfixable_row_indices
  )
  out$unfixable <- as.integer(seq_len(n) %in% assessment$unfixable_row_indices)
  out <- out[, all_cols]

  if (isTRUE(verbose)) {
    n_fix_rows <- sum(out$fixable & !out$unfixable)
    n_unfix_rows <- sum(out$unfixable)
    if (n_fix_rows == 0 && n_unfix_rows == 0) {
      message("  All clean.")
    } else {
      message("  Fixable:   ", n_fix_rows)
      message("  Unfixable: ", n_unfix_rows)
    }
    fix_counts <- vapply(fixable_issue_cols, \(col) sum(out[[col]]), integer(1))
    fix_counts <- fix_counts[fix_counts > 0]
    unfix_counts <- vapply(
      unfixable_issue_cols,
      \(col) sum(out[[col]]),
      integer(1)
    )
    unfix_counts <- unfix_counts[unfix_counts > 0]
    if (length(fix_counts) > 0) {
      message(
        "  Fixable breakdown:   ",
        paste(names(fix_counts), fix_counts, sep = ": ", collapse = ", ")
      )
    }
    if (length(unfix_counts) > 0) {
      message(
        "  Unfixable breakdown: ",
        paste(names(unfix_counts), unfix_counts, sep = ": ", collapse = ", ")
      )
    }
  }

  if (!is.null(id_col_name)) {
    out[[id_col_name]] <- id_values
    out <- dplyr::relocate(out, dplyr::all_of(id_col_name), .before = 1)
  }

  class(out) <- c("karyo_check", class(out))
  attr(out, ".kp_assessment") <- assessment
  attr(out, ".kp_raw") <- karyotypes
  attr(out, ".kp_karyotype_col") <- karyotype_col_name
  attr(out, ".kp_id_col") <- id_col_name
  attr(out, ".kp_id_values") <- id_values

  out
}
