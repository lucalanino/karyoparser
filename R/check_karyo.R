#' Check Karyotype Strings for Issues
#'
#' Checks raw karyotype strings for formatting artifacts and structural errors.
#' Returns a wide-format tibble with one row per input string: each possible
#' issue type is a column (`0`/`1`), plus summary `fixable` and `unfixable`
#' columns. When `verbose = TRUE`, prints a count summary and per-issue-type
#' breakdown.
#'
#' Fixable issues (resolvable by `parse_karyo()` under
#' `on_issues = "preprocess"`):
#' dirty markers (`unicode_notation`, `embedded_newline`, `html_entities`,
#' `leading_dot`, `fish_notation`, `trailing_narrative`, `midstring_linewrap`,
#' `missing_sex_comma`, `mar_space`) and `chimeric_separator`.
#'
#' Unfixable issues (always returned as NA by `parse_karyo()`):
#' `zero_host_chimera`, `empty`, `no_chromosome_count`, `updated_iscn`,
#' `unbalanced_parentheses`, `unbalanced_brackets`, `no_sex_complement`,
#' `invalid_idem`, `unparseable_bracket`.
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
#'   common names; in an interactive session the detected column is shown with a
#'   preview and confirmed before use. Ignored when `karyotypes` is a character
#'   vector.
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

# Run the full fix-and-classify pipeline; shared by check, preprocess, and parse.
.assess_karyotypes <- function(x) {
  DIRTY_TYPES <- setdiff(.fixable_issue_types, "chimeric_separator")

  dirty_issues <- flag_unpreprocessed(x)
  dirty_row_indices <- unique(
    dirty_issues$row_index[dirty_issues$issue_type %in% DIRTY_TYPES]
  )

  partially_fixed <- .apply_dirty_fixes(x)

  # Mirrors validate_karyotypes(): must start with a digit and contain '//' after normalization.
  norm_partial <- normalize_iscn(partially_fixed)
  chimeric_row_indices <- which(
    !is.na(norm_partial) &
      stringr::str_detect(norm_partial, "^\\d") &
      stringr::str_detect(norm_partial, "//")
  )
  chimeric_issues <- if (length(chimeric_row_indices) > 0) {
    tibble::tibble(
      row_index = chimeric_row_indices,
      karyotype = truncate_str(as.character(partially_fixed[
        chimeric_row_indices
      ])),
      issue_type = "chimeric_separator",
      issue_detail = "Contains '//' chimeric separator (independent cell populations)"
    )
  } else {
    empty_issues_tibble()
  }

  fully_fixed <- partially_fixed
  if (length(chimeric_row_indices) > 0) {
    fully_fixed[chimeric_row_indices] <- stringr::str_replace(
      fully_fixed[chimeric_row_indices],
      "//.*$",
      ""
    )
  }

  norm_fully_fixed <- normalize_iscn(fully_fixed)
  struct_issues <- validate_karyotypes(norm_fully_fixed)

  reported_issues <- dplyr::bind_rows(
    dirty_issues,
    chimeric_issues,
    struct_issues
  ) |>
    dplyr::arrange(row_index)

  unfixable_types <- setdiff(.all_issue_types, .fixable_issue_types)
  unfixable_row_indices <- unique(
    reported_issues$row_index[reported_issues$issue_type %in% unfixable_types]
  )

  processed <- norm_fully_fixed
  processed[unfixable_row_indices] <- NA_character_

  list(
    processed = processed,
    reported_issues = reported_issues,
    unfixable_row_indices = unfixable_row_indices,
    dirty_row_indices = dirty_row_indices,
    chimeric_row_indices = chimeric_row_indices
  )
}

# Single source of truth for the dirty-fix loop.
.apply_dirty_fixes <- function(x) {
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "[\n\r\t]+", " ")
  for (nm in names(.dirty_patterns)) {
    dp <- .dirty_patterns[[nm]]
    for (fx in dp$fix) {
      x <- stringr::str_replace_all(x, fx$pattern, fx$replacement)
    }
  }
  stringr::str_trim(x)
}

# Check normalized karyotype strings for structural errors; returns long issues tibble.
validate_karyotypes <- function(karyotypes) {
  issue_list <- vector("list", length(karyotypes))
  n_issues <- 0L

  add_issue <- function(idx, karyo, type, detail) {
    n_issues <<- n_issues + 1L
    issue_list[[n_issues]] <<- tibble::tibble(
      row_index = idx,
      karyotype = truncate_str(as.character(karyo)),
      issue_type = type,
      issue_detail = detail
    )
  }

  for (i in seq_along(karyotypes)) {
    k <- karyotypes[i]

    if (is.na(k) || k == "") {
      add_issue(i, if (is.na(k)) "NA" else "", "empty", "NA or empty string")
      next
    }

    if (!stringr::str_detect(k, "^\\d")) {
      add_issue(
        i,
        k,
        "no_chromosome_count",
        "String doesn't start with chromosome count"
      )
      next
    }

    if (stringr::str_detect(k, "//")) {
      add_issue(
        i,
        k,
        "chimeric_separator",
        "Contains '//' chimeric separator (independent cell populations)"
      )
      next
    }

    if (stringr::str_detect(k, "(?i)Updated ISCN")) {
      add_issue(
        i,
        k,
        "updated_iscn",
        "Contains 'Updated ISCN' correction marker — original karyotype may be superseded"
      )
      next
    }

    open_parens <- stringr::str_count(k, "\\(")
    close_parens <- stringr::str_count(k, "\\)")
    if (open_parens != close_parens) {
      add_issue(
        i,
        k,
        "unbalanced_parentheses",
        sprintf(
          "Mismatched parentheses: %d open, %d close",
          open_parens,
          close_parens
        )
      )
      next
    }

    open_brackets <- stringr::str_count(k, "\\[")
    close_brackets <- stringr::str_count(k, "\\]")
    if (open_brackets != close_brackets) {
      add_issue(
        i,
        k,
        "unbalanced_brackets",
        sprintf(
          "Mismatched brackets: %d open, %d close",
          open_brackets,
          close_brackets
        )
      )
      next
    }

    first_clone <- stringr::str_split(k, "/")[[1]][1]
    first_clone_clean <- stringr::str_replace_all(
      first_clone,
      "\\[[^\\]]+\\]",
      ""
    )
    tokens <- stringr::str_split(first_clone_clean, ",")[[1]]
    has_sex <- any(tokens %in% .sex_complements)
    if (!has_sex && length(tokens) > 1) {
      add_issue(
        i,
        k,
        "no_sex_complement",
        "No sex chromosome complement found in first clone"
      )
    }

    if (length(tokens) > 1 && any(tolower(tokens) == "idem")) {
      add_issue(
        i,
        k,
        "invalid_idem",
        "idem found in clone 1 (nothing to inherit from)"
      )
    }

    brackets <- stringr::str_extract_all(k, "\\[[^\\]]+\\]")[[1]]
    for (br in brackets) {
      content <- stringr::str_replace_all(br, "\\[|\\]", "")
      is_valid_bracket <- stringr::str_detect(content, "^(cp)?\\d+(~\\d+)?$")
      if (!is_valid_bracket && content != "") {
        add_issue(
          i,
          k,
          "unparseable_bracket",
          sprintf("Cannot parse '%s' as metaphase count", br)
        )
      }
    }
  }

  if (n_issues == 0L) {
    return(empty_issues_tibble())
  }
  dplyr::bind_rows(issue_list[seq_len(n_issues)])
}
