#' @include karyoparser-package.R
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
  # Non-overlapping character-level substitutions, so all entries are
  # collapsed into one named-vector fix step (a single str_replace_all() call
  # matching every entry at once) rather than one sequential pass per entry.
  unicode_notation = list(
    detect = "[\u00A0\u2007\u202F\u200B\u2212\u2012\u2013\u2014\uFE63\uFF0D\uFF0B]",
    use_trimmed = FALSE,
    fix = list(c(
      "[\u00A0\u2007\u202F]" = " ",
      "\u200B" = "",
      "[\u2212\u2012\u2013\u2014\uFE63\uFF0D]" = "-",
      "\uFF0B" = "+"
    )),
    detail = "Contains unicode spaces (NBSP), zero-width spaces, dashes (em/en-dash), or fullwidth characters"
  ),
  non_ascii_homoglyph = list(
    detect = "[\u03A7\u03A5]",
    use_trimmed = FALSE,
    fix = list(c("\u03A7" = "X", "\u03A5" = "Y")),
    detail = "Greek letter homoglyph (Chi/Upsilon) standing in for Latin X/Y in sex chromosome complement"
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
    fix = list(c("&lt;" = "<", "&gt;" = ">", "&amp;" = "&")),
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
  # Lookahead includes +/- (not just ,[/ or end) so a count glued to a sex
  # complement that's *also* glued to the first aberration (e.g. '46XY+13')
  # still gets its count-sex comma inserted here; missing_sex_comma (below)
  # then closes the remaining gap on the same pass.
  count_sex_separator = list(
    detect = paste0(
      "(?:^|/)\\d+(?:~\\d+)?[.]?(?:",
      .sex_alt,
      ")(?=[,\\[/+-]|$)"
    ),
    use_trimmed = TRUE,
    fix = list(
      list(
        pattern = paste0(
          "(^|/)(\\d+(?:~\\d+)?)[.]?(",
          .sex_alt,
          ")(?=[,\\[/+-]|$)"
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
    detect = c(
      paste0(",(", .sex_alt, ")\\s+(?=[a-z(+])"),
      # No separator at all, e.g. '46,XX+8' or '46,XY-7' -- covers every
      # clone in the string (not just the first), since the pattern only
      # anchors on a preceding comma, not string start.
      paste0(",(", .sex_alt, ")(?=[+-])"),
      # Glued directly to a letter/paren-led aberration token with no
      # separator at all, e.g. '46,XYdel(5q)'. Detect-only: unlike the two
      # cases above, there's no fix for this below (see TODO.md) -- '[a-z(]'
      # right after a sex complement is too broad a trigger surface to safely
      # auto-insert a comma without real-corpus validation. The negative
      # lookahead excludes the constitutional-sex-complement 'c' suffix (e.g.
      # '47,XXYc[20]', optionally '?'-suffixed) so this doesn't collide with
      # `constitutional_sex_complement`.
      paste0(",(", .sex_alt, ")(?!c\\??(?:[,/\\[]|$))(?=[a-z(])")
    ),
    use_trimmed = FALSE,
    fix = list(
      list(
        pattern = paste0("(,(", .sex_alt, "))\\s+(?=[a-z(+])"),
        replacement = "\\1,"
      ),
      list(
        pattern = paste0("(,(", .sex_alt, "))(?=[+-])"),
        replacement = "\\1,"
      )
    ),
    detail = "Missing comma between sex chromosome complement and first aberration (e.g. '46,XX der(...)' or '46,XY+8' should be '46,XX,der(...)' / '46,XY,+8'; a letter-glued form like '46,XYdel(5q)' is detected but not auto-fixed)"
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
  ),
  non_clonal_sca = list(
    detect = "(?i)\\bncSCA\\b",
    use_trimmed = FALSE,
    fix = list(
      list(pattern = "(?i)\\(ncSCA\\)", replacement = ""),
      list(pattern = "(?i)ncSCA(?:\\[[^\\]]*\\])?", replacement = "")
    ),
    detail = "Non-clonal single-cell abnormality ('ncSCA' token); stripped, remaining clone parsed"
  ),
  # Last so it only ever sees what survives every earlier fix (including
  # trailing_narrative's stripping) -- catches leftover non-ASCII with no
  # known automatic fix, without false-firing on discarded narrative text.
  stray_non_ascii = list(
    detect = "[^\x01-\x7F]",
    use_trimmed = FALSE,
    fix = list(),
    detail = "Contains non-ASCII character(s) not recognized by any automatic fix; manual review needed"
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
  "single_token",
  "constitutional_sex_complement",
  "mosaic_karyotype",
  "invalid_idem",
  "unparseable_bracket"
)

# Single ordered walk over .dirty_patterns, replacing what used to be two
# independent passes (one for detection, one for fixing) over the raw string.
# Threads a "state" vector through the patterns in list order: each pattern's
# detect runs against the state as left by every earlier pattern's fix (not
# the raw original), then that pattern's own fix is applied before moving to
# the next one. This mirrors the dependency the fixes already have on each
# other (e.g. fish_notation must claim a trailing FISH clause before
# trailing_narrative's blunter truncation rule would misclassify it as
# narrative) so detection and fixing finally agree on what each row's issues
# are.
.run_dirty_patterns <- function(x) {
  orig_display <- truncate_str(as.character(x))
  na_mask <- is.na(x)

  state <- stringr::str_trim(x)
  issue_list <- vector("list", length(.dirty_patterns))

  for (i in seq_along(.dirty_patterns)) {
    nm <- names(.dirty_patterns)[i]
    dp <- .dirty_patterns[[i]]

    target <- if (isTRUE(dp$use_trimmed)) stringr::str_trim(state) else state
    detected <- Reduce(
      `|`,
      lapply(dp$detect, function(re) {
        hit <- stringr::str_detect(target, re)
        hit[is.na(hit)] <- FALSE
        hit
      })
    )
    detected[na_mask] <- FALSE

    hit_rows <- which(detected)
    if (length(hit_rows) > 0) {
      issue_list[[i]] <- tibble::tibble(
        row_index = hit_rows,
        karyotype = orig_display[hit_rows],
        issue_type = nm,
        issue_detail = dp$detail
      )
    }

    for (fx in dp$fix) {
      # Gated on this pattern's own `detected`: a permissive fix regex (e.g.
      # trailing_narrative's "everything after the last bracket" truncation)
      # must not touch rows it never reported as having this issue, or it
      # silently destroys content -- including a later clone's structural
      # errors -- before a downstream pattern gets a chance to flag it.
      #
      # A step is either list(pattern=, replacement=) -- one sequential pass,
      # for fixes whose patterns depend on an earlier step's output -- or a
      # named character vector (pattern = replacement) applying several
      # non-overlapping character-level substitutions in a single pass.
      fixed_state <- if (is.list(fx)) {
        stringr::str_replace_all(state, fx$pattern, fx$replacement)
      } else {
        stringr::str_replace_all(state, fx)
      }
      state <- ifelse(detected, fixed_state, state)
    }
  }

  issues <- dplyr::bind_rows(Filter(Negate(is.null), issue_list))
  if (nrow(issues) == 0L) {
    issues <- empty_issues_tibble()
  }

  list(issues = issues, fixed = stringr::str_trim(state))
}

# Aggregates a flat per-token/per-item logical vector back to one value per
# row (TRUE if any item in that row's group is TRUE). Assumes every group in
# seq_len(n_groups) has at least one member (true of str_split() output,
# which always returns >= 1 element per input string).
.agg_any_by_group <- function(flat_logical, flat_group, n_groups) {
  as.logical(tapply(
    flat_logical,
    factor(flat_group, levels = seq_len(n_groups)),
    any
  ))
}

validate_karyotypes <- function(karyotypes) {
  n <- length(karyotypes)
  disp <- truncate_str(as.character(karyotypes))

  # Constitutional sex complement, e.g. '47,XXYc'; out of scope, flagged unfixable
  # rather than mislabelled no_sex_complement.
  const_sex_re <- paste0("^(?:", .sex_alt, ")c\\??$")

  # An aberrant sex chromosome (e.g. '-Y', '+X', or X/Y inside a der/t() paren list)
  # still accounts for sex, so it must not trip no_sex_complement.
  sex_aberr_re <- "(?:^[+-](?:X|Y)(?:x\\d+)?$)|[(;](?:X|Y)[);]"

  issue_list <- list()
  add_batch <- function(idx, karyo, type, detail) {
    if (length(idx) == 0L) {
      return(invisible())
    }
    issue_list[[length(issue_list) + 1L]] <<- tibble::tibble(
      row_index = idx,
      karyotype = karyo,
      issue_type = type,
      issue_detail = detail
    )
  }

  # Priority chain: each check below only claims rows still `remaining`, so a
  # row is reported under (at most) the first category that matches it, same
  # as the `next`-chained early exits this replaces. `remaining & <NA-capable
  # expr>` is always FALSE (never NA) for rows already excluded, since `&`
  # short-circuits to FALSE on a FALSE left side.
  remaining <- rep(TRUE, n)

  is_na <- is.na(karyotypes)
  hit <- remaining & (is_na | karyotypes == "")
  idx <- which(hit)
  add_batch(idx, ifelse(is_na[idx], "NA", ""), "empty", "NA or empty string")
  remaining <- remaining & !hit

  # Checked before the count/sex checks so these out-of-scope constructs
  # aren't mislabelled no_chromosome_count / no_sex_complement.
  hit <- remaining & stringr::str_detect(karyotypes, "(?i)^mos\\b")
  idx <- which(hit)
  add_batch(
    idx,
    disp[idx],
    "mosaic_karyotype",
    "Mosaic 'mos' prefix; mosaic karyotypes are out of scope"
  )
  remaining <- remaining & !hit

  hit <- remaining & !stringr::str_detect(karyotypes, "^\\d")
  idx <- which(hit)
  add_batch(
    idx,
    disp[idx],
    "no_chromosome_count",
    "String doesn't start with chromosome count"
  )
  remaining <- remaining & !hit

  hit <- remaining & stringr::str_detect(karyotypes, "//")
  idx <- which(hit)
  add_batch(
    idx,
    disp[idx],
    "chimeric_separator",
    "Contains '//' chimeric separator (independent cell populations)"
  )
  remaining <- remaining & !hit

  hit <- remaining & stringr::str_detect(karyotypes, "(?i)Updated ISCN")
  idx <- which(hit)
  add_batch(
    idx,
    disp[idx],
    "updated_iscn",
    "Contains 'Updated ISCN' correction marker \u2014 original karyotype may be superseded"
  )
  remaining <- remaining & !hit

  open_parens <- stringr::str_count(karyotypes, "\\(")
  close_parens <- stringr::str_count(karyotypes, "\\)")
  hit <- remaining & (open_parens != close_parens)
  idx <- which(hit)
  add_batch(
    idx,
    disp[idx],
    "unbalanced_parentheses",
    sprintf(
      "Mismatched parentheses: %d open, %d close",
      open_parens[idx],
      close_parens[idx]
    )
  )
  remaining <- remaining & !hit

  open_brackets <- stringr::str_count(karyotypes, "\\[")
  close_brackets <- stringr::str_count(karyotypes, "\\]")
  hit <- remaining & (open_brackets != close_brackets)
  idx <- which(hit)
  add_batch(
    idx,
    disp[idx],
    "unbalanced_brackets",
    sprintf(
      "Mismatched brackets: %d open, %d close",
      open_brackets[idx],
      close_brackets[idx]
    )
  )
  remaining <- remaining & !hit

  # From here on, only rows that survived every check above are in play, and
  # tokens/brackets are ragged (variable count per row) -- flatten to (row,
  # item) pairs across the whole remaining subset so each regex runs once
  # over every item in the dataset, instead of once per item per row.
  rem_idx <- which(remaining)
  if (length(rem_idx) > 0L) {
    kk <- karyotypes[rem_idx]

    first_clone <- vapply(stringr::str_split(kk, "/"), `[`, character(1), 1)
    first_clone_clean <- stringr::str_replace_all(
      first_clone,
      "\\[[^\\]]+\\]",
      ""
    )
    tokens_list <- stringr::str_split(first_clone_clean, ",")
    n_tokens <- lengths(tokens_list)
    flat_tokens <- unlist(tokens_list)
    flat_row <- rep(seq_along(tokens_list), n_tokens)

    has_sex <- .agg_any_by_group(
      flat_tokens %in% .sex_complements,
      flat_row,
      length(tokens_list)
    )
    has_const_sex <- .agg_any_by_group(
      stringr::str_detect(flat_tokens, const_sex_re),
      flat_row,
      length(tokens_list)
    )
    has_sex_aberr <- .agg_any_by_group(
      stringr::str_detect(flat_tokens, sex_aberr_re),
      flat_row,
      length(tokens_list)
    )
    has_idem <- .agg_any_by_group(
      tolower(flat_tokens) == "idem",
      flat_row,
      length(tokens_list)
    )

    sub_remaining <- rep(TRUE, length(rem_idx))

    hit <- sub_remaining & has_const_sex
    idx <- rem_idx[hit]
    add_batch(
      idx,
      disp[idx],
      "constitutional_sex_complement",
      "Constitutional sex complement (e.g. '47,XXYc'); constitutional abnormalities are out of scope"
    )
    sub_remaining <- sub_remaining & !hit

    hit <- sub_remaining & (n_tokens == 1L)
    idx <- rem_idx[hit]
    add_batch(
      idx,
      disp[idx],
      "single_token",
      "Single token with no comma: ambiguous between a bare chromosome count and an abnormality-only entry with a stripped +/- sign (e.g. Excel autocorrect turning '+8' into '8')"
    )
    sub_remaining <- sub_remaining & !hit

    hit <- sub_remaining & !has_sex & !has_sex_aberr & (n_tokens > 1L)
    idx <- rem_idx[hit]
    add_batch(
      idx,
      disp[idx],
      "no_sex_complement",
      "No sex chromosome complement found in first clone"
    )

    idx <- rem_idx[n_tokens > 1L & has_idem]
    add_batch(
      idx,
      disp[idx],
      "invalid_idem",
      "idem found in clone 1 (nothing to inherit from)"
    )

    brackets_list <- stringr::str_extract_all(kk, "\\[[^\\]]+\\]")
    br_lengths <- lengths(brackets_list)
    if (sum(br_lengths) > 0L) {
      flat_brackets <- unlist(brackets_list)
      flat_br_row <- rem_idx[rep(seq_along(rem_idx), br_lengths)]
      content <- stringr::str_replace_all(flat_brackets, "\\[|\\]", "")
      invalid <- !stringr::str_detect(content, "^(cp)?\\d+(~\\d+)?$") &
        content != ""
      add_batch(
        flat_br_row[invalid],
        disp[flat_br_row[invalid]],
        "unparseable_bracket",
        sprintf("Cannot parse '%s' as metaphase count", flat_brackets[invalid])
      )
    }
  }

  if (length(issue_list) == 0L) {
    return(empty_issues_tibble())
  }
  dplyr::bind_rows(issue_list)
}

# Shared fix-and-classify pipeline for check/preprocess/parse.
.assess_karyotypes <- function(x, on_chimeric = c("default", "host", "donor")) {
  on_chimeric <- match.arg(on_chimeric)
  n <- length(x)
  DIRTY_TYPES <- setdiff(
    .fixable_issue_types,
    c("chimeric_separator", "zero_host_chimera")
  )

  dp_walk <- .run_dirty_patterns(x)
  dirty_issues <- dp_walk$issues
  dirty_row_indices <- unique(
    dirty_issues$row_index[dirty_issues$issue_type %in% DIRTY_TYPES]
  )
  zhc_row_indices <- unique(
    dirty_issues$row_index[dirty_issues$issue_type == "zero_host_chimera"]
  )

  # zero_host_chimera's fix already strips the leading './/', reducing it to a donor.
  partially_fixed <- dp_walk$fixed
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
