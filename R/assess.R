#' @include parse_karyo-package.R
NULL

truncate_str <- function(s, max_len = 40) {
  ifelse(nchar(s) > max_len, paste0(substr(s, 1, max_len - 3), "..."), s)
}

empty_issues_tibble <- function() {
  tibble::tibble(
    row_index = integer(),
    karyotype = character(),
    issue_type = character(),
    issue_detail = character()
  )
}

# Longest-first so the regex engine can't match a shorter prefix before a longer one.
.sex_alt <- paste(
  .sex_complements[order(-nchar(.sex_complements))],
  collapse = "|"
)

# detect/fix/detail entries; empty fix list means detected-only (unfixable).
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
    detect = "^[.\\s]*[.]\\s*(?=\\d)",
    use_trimmed = TRUE,
    fix = list(
      list(pattern = "^[.\\s]*[.]\\s*(?=\\d)", replacement = "")
    ),
    detail = "String starts with dot(s) (and any stray whitespace) before chromosome count"
  ),
  # Runs after leading_dot; repairs a missing/dotted count-sex separator, e.g.
  # '46XY,...' or '45.XY,...' -> '46,XY,...'.
  count_sex_separator = list(
    detect = paste0(
      "(?:^|/)\\d+(?:~\\d+)?[.]?(?:",
      .sex_alt,
      ")(?=[,\\[/]|$)"
    ),
    use_trimmed = TRUE,
    fix = list(
      list(
        pattern = paste0(
          "(^|/)(\\d+(?:~\\d+)?)[.]?(",
          .sex_alt,
          ")(?=[,\\[/]|$)"
        ),
        replacement = "\\1\\2,\\3"
      )
    ),
    detail = "Missing or dotted separator between chromosome count and sex complement (e.g. '46XY' or '45.XY' should be '46,XY')"
  ),
  fish_notation = list(
    detect = "[. ]+(?:nuc )?ish\\b",
    use_trimmed = FALSE,
    fix = list(
      # 'nuc ish' brackets are FISH cell counts, not the metaphase count; strip
      # the whole clause to end of string (or a '//' chimeric boundary).
      list(
        pattern = "[. ]+nuc ish\\b.*?(?=//|$)",
        replacement = ""
      ),
      # Bare '.ish'/'ish' may sit mid-clone; stop before the metaphase count bracket.
      list(
        pattern = "[. ]+ish\\b.*?(?=\\[(?:cp)?\\d+(?:[~-]\\d+)?\\]|$)",
        replacement = ""
      )
    ),
    detail = "FISH / nuc ish annotation (suffix or mid-clone before metaphase count bracket)"
  ),
  trailing_narrative = list(
    detect = c(
      "\\]\\s*\\.",
      "\\]\\s+[A-Z][a-z]",
      "\\s+\\.\\s*[A-Z]",
      "\\)\\s+[A-Z][a-z]",
      # Narrative after a count+sex normal karyotype with no bracket/paren.
      paste0("^\\d+(?:~\\d+)?,(?:", .sex_alt, ")\\s+[A-Z][a-z]")
    ),
    use_trimmed = FALSE,
    fix = list(
      # Strip everything after the last metaphase-count bracket ([n]/[cpN]/[n~m]).
      list(
        pattern = "^(.*\\[(?:cp)?\\d+(?:[~-]\\d+)?\\]).*$",
        replacement = "\\1"
      ),
      # Collapse ] ./ separator artifact left between clones.
      list(
        pattern = "\\]\\s*\\./",
        replacement = "]/"
      ),
      # Strip narrative after last ')' when no metaphase bracket is present.
      list(
        pattern = "^(.*\\))\\s+[A-Z][a-z].*$",
        replacement = "\\1"
      ),
      # Fallback for strings with no metaphase bracket or parens.
      list(
        pattern = "\\s+\\.\\s*[A-Z].*$",
        replacement = ""
      ),
      # Strip narrative after a count+sex karyotype with no bracket/paren to anchor
      # on; the lookahead keeps an "Updated ISCN" marker unfixable rather than
      # silently truncating it to a clean count+sex.
      list(
        pattern = paste0(
          "^(?!.*(?i:Updated ISCN))(\\d+(?:~\\d+)?,(?:",
          .sex_alt,
          "))\\s+[A-Z][a-z].*$"
        ),
        replacement = "\\1"
      )
    ),
    detail = "Trailing narrative text after last metaphase-count bracket, closing parenthesis, or count+sex normal karyotype"
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
  # Detected before normalize_iscn() strips leading punctuation.
  zero_host_chimera = list(
    detect = "^[.]*//",
    use_trimmed = TRUE,
    fix = list(
      list(pattern = "^[.]*//", replacement = "")
    ),
    detail = "String starts with './/' or '//': donor-only chimera with no host metaphases"
  )
)

# Defined after .dirty_patterns so load order is guaranteed.
.fixable_issue_types <- c(
  names(Filter(\(p) length(p$fix) > 0, .dirty_patterns)),
  "chimeric_separator"
)

.all_issue_types <- c(
  names(.dirty_patterns),
  "empty",
  "no_chromosome_count",
  "chimeric_separator",
  "multiple_chimeric_separator",
  "updated_iscn",
  "unbalanced_parentheses",
  "unbalanced_brackets",
  "no_sex_complement",
  "constitutional_sex_complement",
  "mosaic_karyotype",
  "non_clonal_sca",
  "invalid_idem",
  "unparseable_bracket"
)

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

  # Constitutional sex complement, e.g. '47,XXYc'; out of scope, flagged unfixable
  # rather than mislabelled no_sex_complement.
  const_sex_re <- paste0("^(?:", .sex_alt, ")c\\??$")

  # An aberrant sex chromosome (e.g. '-Y', '+X', or X/Y inside a der/t() paren list)
  # still accounts for sex, so it must not trip no_sex_complement.
  sex_aberr_re <- "(?:^[+-](?:X|Y)(?:x\\d+)?$)|[(;](?:X|Y)[);]"

  for (i in seq_along(karyotypes)) {
    k <- karyotypes[i]

    if (is.na(k) || k == "") {
      add_issue(i, if (is.na(k)) "NA" else "", "empty", "NA or empty string")
      next
    }

    # Checked before the count/sex checks so these out-of-scope constructs
    # aren't mislabelled no_chromosome_count / no_sex_complement.
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

# Shared fix-and-classify pipeline for check/preprocess/parse.
.assess_karyotypes <- function(x, on_chimeric = c("default", "host", "donor")) {
  on_chimeric <- match.arg(on_chimeric)
  n <- length(x)
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

  # zero_host_chimera's fix already strips the leading './/', reducing it to a donor.
  partially_fixed <- .apply_dirty_fixes(x)
  norm_partial <- normalize_iscn(partially_fixed)

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

  if (length(multi_row_indices) > 0) {
    selected[multi_row_indices] <- NA_character_
  }

  struct_issues <- validate_karyotypes(selected)
  # Rows blanked by chimeric policy must not surface as structural 'empty'.
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

  # Zero-host-under-"host" NAs are policy-induced, not a data defect, so the row
  # stays classified fixable (zero_host_chimera), not unfixable.
  unfixable_types <- setdiff(.all_issue_types, .fixable_issue_types)
  unfixable_row_indices <- unique(
    reported_issues$row_index[reported_issues$issue_type %in% unfixable_types]
  )

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
