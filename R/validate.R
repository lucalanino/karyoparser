#' Check Karyotype Strings for Issues
#'
#' Runs both dirty-marker detection and structural validation on raw karyotype
#' strings. Returns a tibble of all issues found, suitable for inspection before
#' parsing with `parse_karyo()`.
#'
#' @param x Character vector of karyotype strings.
#' @return A tibble with columns: `row_index`, `karyotype` (truncated to 40
#'   chars), `issue_type`, `issue_detail`. Each issue generates one row;
#'   multiple issues on the same input row appear as separate rows.
#' @export
check_karyo <- function(x) {
  if (length(x) == 0) {
    return(tibble::tibble(
      row_index = integer(),
      karyotype = character(),
      issue_type = character(),
      issue_detail = character()
    ))
  }
  dirty_issues <- flag_unpreprocessed(x)
  normalized <- normalize_iscn(x)
  struct_issues <- validate_karyotypes(normalized)
  dplyr::bind_rows(dirty_issues, struct_issues) |>
    dplyr::arrange(row_index)
}

#' Validate Karyotype Strings
#'
#' Checks karyotype strings for common structural issues. Called internally by
#' `check_karyo()` after normalization.
#'
#' @param karyotypes Character vector of karyotype strings (already preprocessed)
#' @return A tibble with columns: row_index, karyotype (truncated), issue_type, issue_detail
#' @keywords internal
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

    # Check 1: NA or empty
    if (is.na(k) || k == "") {
      add_issue(i, ifelse(is.na(k), "NA", ""), "empty", "NA or empty string")
      next
    }

    # Check 2: No leading chromosome count
    if (!stringr::str_detect(k, "^\\d")) {
      add_issue(
        i,
        k,
        "no_chromosome_count",
        "String doesn't start with chromosome count"
      )
      next
    }

    # Check 2b: Chimeric separator (//)
    if (stringr::str_detect(k, "//")) {
      add_issue(
        i,
        k,
        "chimeric_separator",
        "Contains '//' chimeric separator (independent cell populations)"
      )
      next
    }

    # Check 2c: Updated ISCN correction marker
    if (stringr::str_detect(k, "(?i)Updated ISCN")) {
      add_issue(
        i,
        k,
        "updated_iscn",
        "Contains 'Updated ISCN' correction marker \u2014 original karyotype may be superseded"
      )
      next
    }

    # Check 3: Unbalanced parentheses
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

    # Check 4: Unbalanced brackets
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

    # Check 5: No sex complement in first clone
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

    # Check 6: idem in clone 1
    if (length(tokens) > 1 && any(tolower(tokens) == "idem")) {
      add_issue(
        i,
        k,
        "invalid_idem",
        "idem found in clone 1 (nothing to inherit from)"
      )
    }

    # Check 7: Unparseable bracket content
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
    return(tibble::tibble(
      row_index = integer(),
      karyotype = character(),
      issue_type = character(),
      issue_detail = character()
    ))
  }
  dplyr::bind_rows(issue_list[seq_len(n_issues)])
}

#' Print Issue Report
#' @param issues Tibble returned by check_karyo()
#' @param total Total number of karyotypes checked
#' @keywords internal
print_issue_report <- function(issues, total) {
  n_issues <- length(unique(issues$row_index))
  n_valid <- total - n_issues

  cat("\n")
  cat("Karyotype Issue Report\n")
  cat("======================\n")
  cat(sprintf("Checked: %d karyotypes\n", total))
  cat(sprintf("Clean:   %d (%.1f%%)\n", n_valid, 100 * n_valid / total))
  cat(sprintf("Issues:  %d (%.1f%%)\n", n_issues, 100 * n_issues / total))

  if (nrow(issues) > 0) {
    cat("\nIssues found:\n")

    # Group by issue type for cleaner output
    by_type <- issues |>
      dplyr::group_by(issue_type) |>
      dplyr::summarise(
        count = dplyr::n(),
        rows = paste(unique(row_index), collapse = ", "),
        .groups = "drop"
      )

    for (j in seq_len(nrow(by_type))) {
      cat(sprintf(
        "  - %s (%d): rows %s\n",
        by_type$issue_type[j],
        by_type$count[j],
        by_type$rows[j]
      ))
    }

    cat("\nDetailed issues:\n")
    print(issues, n = min(20, nrow(issues)))
    if (nrow(issues) > 20) {
      cat(sprintf("  ... and %d more issues\n", nrow(issues) - 20))
    }
  }
  cat("\n")
}
