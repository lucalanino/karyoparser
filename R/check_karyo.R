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
#' `missing_sex_comma`, `mar_space`), `chimeric_separator`, and
#' `zero_host_chimera` (the donor clone is parsed; see `on_chimeric` in
#' `parse_karyo()`).
#'
#' Unfixable issues (always returned as NA by `parse_karyo()`):
#' `multiple_chimeric_separator`, `empty`, `no_chromosome_count`,
#' `updated_iscn`, `unbalanced_parentheses`, `unbalanced_brackets`,
#' `no_sex_complement`, `constitutional_sex_complement`, `mosaic_karyotype`,
#' `non_clonal_sca`, `invalid_idem`, `unparseable_bracket`.
#'
#' Some constructs are deliberately out of scope and flagged unfixable rather
#' than parsed: constitutional abnormalities (a sex complement with a `c`
#' suffix, e.g. `47,XXYc`, flagged `constitutional_sex_complement`), mosaic
#' karyotypes (a leading `mos` prefix, flagged `mosaic_karyotype`), and
#' non-clonal single-cell abnormalities (an `ncSCA` token, flagged
#' `non_clonal_sca`). These take precedence over `no_chromosome_count` and
#' `no_sex_complement` so a present count or complement is not mislabelled.
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

# Run the full fix-and-classify pipeline; shared by check, preprocess, and parse.
# Clone selection for chimeric rows is governed by `on_chimeric`:
#   "default" -- host clone for normal chimeras, donor for zero-host chimeras
#   "host"    -- host clone only; zero-host chimeras become NA (unfixable)
#   "donor"   -- everything after '//' (the donor population)
.assess_karyotypes <- function(x, on_chimeric = c("default", "host", "donor")) {
  on_chimeric <- match.arg(on_chimeric)
  n <- length(x)
  # zero_host_chimera is a chimeric concern, not a plain dirty fix.
  DIRTY_TYPES <- setdiff(
    .fixable_issue_types,
    c("chimeric_separator", "zero_host_chimera")
  )

  dirty_issues <- flag_unpreprocessed(x)
  dirty_row_indices <- unique(
    dirty_issues$row_index[dirty_issues$issue_type %in% DIRTY_TYPES]
  )
  zhc_row_indices <- unique(
    dirty_issues$row_index[dirty_issues$issue_type == "zero_host_chimera"]
  )

  # Dirty fixes include the zero_host_chimera fix, which strips a leading './/'
  # prefix -- so zero-host rows arrive here already reduced to their donor.
  partially_fixed <- .apply_dirty_fixes(x)
  norm_partial <- normalize_iscn(partially_fixed)

  # Classify chimeric rows by counting '//' separators in the normalized string.
  has_count <- !is.na(norm_partial) & stringr::str_detect(norm_partial, "^\\d")
  n_sep <- ifelse(
    is.na(norm_partial),
    0L,
    stringr::str_count(norm_partial, stringr::fixed("//"))
  )
  multi_row_indices <- setdiff(which(has_count & n_sep >= 2L), zhc_row_indices)
  chimeric_row_indices <- setdiff(
    which(has_count & n_sep == 1L),
    zhc_row_indices
  )

  # ---- Clone selection per on_chimeric --------------------------------------
  selected <- norm_partial
  chimeric_clone <- rep(NA_character_, n)

  if (length(chimeric_row_indices) > 0) {
    if (on_chimeric == "donor") {
      selected[chimeric_row_indices] <- stringr::str_replace(
        norm_partial[chimeric_row_indices],
        "^.*?//",
        ""
      )
      chimeric_clone[chimeric_row_indices] <- "donor"
    } else {
      selected[chimeric_row_indices] <- stringr::str_replace(
        norm_partial[chimeric_row_indices],
        "//.*$",
        ""
      )
      chimeric_clone[chimeric_row_indices] <- "host"
    }
  }

  if (length(zhc_row_indices) > 0) {
    if (on_chimeric == "host") {
      selected[zhc_row_indices] <- NA_character_
    } else {
      chimeric_clone[zhc_row_indices] <- "donor"
    }
  }

  # Two or more '//' separators cannot resolve to a single donor: unfixable.
  if (length(multi_row_indices) > 0) {
    selected[multi_row_indices] <- NA_character_
  }

  # ---- Structural validation of the selected clone --------------------------
  struct_issues <- validate_karyotypes(selected)
  # Rows blanked purely by chimeric policy (donor-only under "host", or 2+ '//')
  # must not be reported as structural 'empty'; they are handled separately.
  policy_na <- c(
    multi_row_indices,
    if (on_chimeric == "host") zhc_row_indices else integer(0)
  )
  if (length(policy_na) > 0 && nrow(struct_issues) > 0) {
    struct_issues <- struct_issues[
      !(struct_issues$row_index %in%
        policy_na &
        struct_issues$issue_type == "empty"),
      ,
      drop = FALSE
    ]
  }

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

  multi_issues <- if (length(multi_row_indices) > 0) {
    tibble::tibble(
      row_index = multi_row_indices,
      karyotype = truncate_str(as.character(partially_fixed[
        multi_row_indices
      ])),
      issue_type = "multiple_chimeric_separator",
      issue_detail = "Contains two or more '//' chimeric separators (cannot resolve a single donor)"
    )
  } else {
    empty_issues_tibble()
  }

  reported_issues <- dplyr::bind_rows(
    dirty_issues,
    chimeric_issues,
    multi_issues,
    struct_issues
  ) |>
    dplyr::arrange(row_index)

  unfixable_types <- setdiff(.all_issue_types, .fixable_issue_types)
  unfixable_row_indices <- unique(c(
    reported_issues$row_index[reported_issues$issue_type %in% unfixable_types],
    # Donor-only chimera under "host" policy yields no parseable clone.
    if (on_chimeric == "host") zhc_row_indices else integer(0)
  ))

  processed <- selected
  processed[unfixable_row_indices] <- NA_character_

  list(
    processed = processed,
    reported_issues = reported_issues,
    unfixable_row_indices = unfixable_row_indices,
    dirty_row_indices = dirty_row_indices,
    chimeric_row_indices = chimeric_row_indices,
    zhc_row_indices = zhc_row_indices,
    multi_row_indices = multi_row_indices,
    chimeric_clone = chimeric_clone,
    on_chimeric = on_chimeric
  )
}

# Single source of truth for the dirty-fix loop.
.apply_dirty_fixes <- function(x) {
  x <- stringr::str_trim(x)
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

  # Sex complement carrying a constitutional 'c' suffix (optionally with a '?'
  # uncertainty marker), e.g. '47,XXYc' or '47,XXXc?'. Constitutional
  # abnormalities are out of scope, so these are flagged unfixable rather than
  # mislabelled as no_sex_complement.
  const_sex_re <- paste0("^(?:", .sex_alt, ")c\\??$")

  # A clone may carry no plain sex complement when its sex chromosomes are
  # themselves aberrant (e.g. '51,add(X)(q26),-Y,...' or '44,-X,t(X;14)...').
  # Such a token still accounts for sex, so it must not trip no_sex_complement.
  # Match a numerical sex gain/loss ('-Y', '+X', '-Xx2') or X/Y appearing as a
  # chromosome inside an aberration's parenthesised list ('(X)', '(X;', ';Y)').
  # The operator/paren context keeps this from matching a bare 'XX(ncSCA)'.
  sex_aberr_re <- "(?:^[+-](?:X|Y)(?:x\\d+)?$)|[(;](?:X|Y)[);]"

  for (i in seq_along(karyotypes)) {
    k <- karyotypes[i]

    if (is.na(k) || k == "") {
      add_issue(i, if (is.na(k)) "NA" else "", "empty", "NA or empty string")
      next
    }

    # Mosaicism ('mos' prefix) and non-clonal single-cell abnormalities
    # ('ncSCA') are out of scope. Flag them specifically -- and before the
    # chromosome-count / sex-complement checks -- so a present count is not
    # mislabelled no_chromosome_count just because of a leading 'mos', and an
    # ncSCA token is not mislabelled no_sex_complement.
    if (stringr::str_detect(k, "(?i)^mos\\b")) {
      add_issue(
        i,
        k,
        "mosaic_karyotype",
        "Mosaic 'mos' prefix; mosaic karyotypes are out of scope"
      )
      next
    }

    if (stringr::str_detect(k, "(?i)ncSCA")) {
      add_issue(
        i,
        k,
        "non_clonal_sca",
        "Non-clonal single-cell abnormalities (ncSCA) are out of scope"
      )
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
        "Contains 'Updated ISCN' correction marker \u2014 original karyotype may be superseded"
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
    has_const_sex <- any(stringr::str_detect(tokens, const_sex_re))
    has_sex_aberr <- any(stringr::str_detect(tokens, sex_aberr_re))
    if (has_const_sex) {
      add_issue(
        i,
        k,
        "constitutional_sex_complement",
        "Constitutional sex complement (e.g. '47,XXYc'); constitutional abnormalities are out of scope"
      )
    } else if (!has_sex && !has_sex_aberr && length(tokens) > 1) {
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
