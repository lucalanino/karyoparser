# ---- Internal pipeline helpers ------------------------------------------------

build_sample_meta <- function(input_df, verbose = FALSE) {
  counts_raw <- input_df |>
    dplyr::mutate(
      bracket = stringr::str_extract_all(original_karyotype, "\\[[^\\]]+\\]")
    ) |>
    tidyr::unnest_longer(bracket, keep_empty = TRUE) |>
    dplyr::mutate(
      content = ifelse(
        is.na(bracket),
        NA_character_,
        stringr::str_replace_all(bracket, "\\[|\\]", "")
      ),
      is_cp = !is.na(content) & stringr::str_detect(content, "^cp\\d+$"),
      is_range_count = !is.na(content) & stringr::str_detect(content, "~"),
      max_count = dplyr::case_when(
        is_cp ~ suppressWarnings(as.numeric(stringr::str_extract(
          content,
          "\\d+"
        ))),
        is_range_count ~ {
          nums <- stringr::str_extract_all(content, "\\d+")
          purrr::map_dbl(
            nums,
            ~ if (length(.x) >= 1) {
              suppressWarnings(max(as.numeric(.x)))
            } else {
              NA_real_
            }
          )
        },
        TRUE ~ suppressWarnings(as.numeric(content))
      )
    )

  unparseable <- counts_raw |>
    dplyr::filter(!is.na(content) & content != "" & is.na(max_count)) |>
    dplyr::distinct(.pk_row_id, content)

  if (nrow(unparseable) > 0 && isTRUE(verbose)) {
    warning(
      sprintf(
        "Could not parse metaphase count from %d bracket(s): %s. These will be treated as 0.",
        nrow(unparseable),
        paste(sprintf("[%s]", unparseable$content), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  counts_tbl <- counts_raw |>
    dplyr::group_by(.pk_row_id) |>
    dplyr::summarise(
      total_metaphases = if (all(is.na(bracket))) {
        NA_integer_
      } else {
        as.integer(sum(tidyr::replace_na(max_count, 0)))
      },
      is_composite = any(is_cp, na.rm = TRUE),
      .groups = "drop"
    )

  input_df |>
    dplyr::mutate(
      head = stringr::str_extract(original_karyotype, "^\\d+(?:~\\d+)?"),
      is_range_karyotype = !is.na(head) & stringr::str_detect(head, "~"),
      cleaned = stringr::str_trim(stringr::str_replace_all(
        original_karyotype,
        "\\[[^\\]]+\\]",
        ""
      )),
      normal_karyotype = ifelse(
        !is_range_karyotype & cleaned %in% c("46,XX", "46,XY"),
        1L,
        0L
      ),
      ploidy_result = purrr::map(original_karyotype, ploidy_category),
      ploidy_category = purrr::map_chr(ploidy_result, "ploidy"),
      chromosome_count = purrr::map_int(ploidy_result, "chromosome_count"),
      mixed_ploidy = as.integer(purrr::map_lgl(ploidy_result, "mixed"))
    ) |>
    dplyr::select(-ploidy_result) |>
    dplyr::left_join(counts_tbl, by = ".pk_row_id")
}

build_clone_tokens <- function(sample_meta) {
  clones_tbl <- sample_meta |>
    dplyr::select(
      .pk_row_id,
      original_karyotype,
      is_composite,
      is_range_karyotype
    ) |>
    dplyr::mutate(
      clone_str = stringr::str_replace_all(
        original_karyotype,
        "\\[[^\\]]+\\]",
        ""
      )
    ) |>
    tidyr::separate_longer_delim(clone_str, delim = "/") |>
    dplyr::group_by(.pk_row_id) |>
    dplyr::mutate(clone_id = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::mutate(clone_str = stringr::str_trim(clone_str)) |>
    tidyr::separate(
      clone_str,
      into = c("chromosome_count", "aberrations"),
      sep = ",",
      fill = "right",
      extra = "merge",
      remove = FALSE
    ) |>
    dplyr::mutate(
      chromosome_count = stringr::str_trim(chromosome_count),
      aberrations = stringr::str_trim(aberrations)
    )

  tokens_tbl <- clones_tbl |>
    dplyr::mutate(aberrations = ifelse(is.na(aberrations), "", aberrations)) |>
    tidyr::separate_longer_delim(aberrations, delim = ",") |>
    dplyr::mutate(aberr_raw = stringr::str_trim(aberrations)) |>
    dplyr::select(
      .pk_row_id,
      clone_id,
      chromosome_count,
      aberr_raw,
      is_range_karyotype,
      is_composite
    ) |>
    dplyr::group_by(.pk_row_id, clone_id) |>
    dplyr::mutate(token_index = dplyr::row_number()) |>
    dplyr::ungroup()

  range_tokens <- sample_meta |>
    dplyr::filter(is_range_karyotype) |>
    dplyr::transmute(
      .pk_row_id,
      clone_id = NA_integer_,
      chromosome_count = NA_character_,
      aberr_raw = "chromosome_count_range",
      token_index = NA_integer_,
      is_range_karyotype = TRUE,
      is_composite = is_composite
    )

  dplyr::bind_rows(tokens_tbl, range_tokens) |>
    dplyr::mutate(aberr_norm = normalize_token(aberr_raw))
}

match_rules <- function(tokens_tbl, rules, rule_flag_names) {
  match_text <- strip_bands(tokens_tbl$aberr_raw)
  n_tokens <- length(match_text)
  n_rules <- nrow(rules)

  # Stage 1: logical matrix (n_tokens x n_rules)
  match_mat <- matrix(FALSE, nrow = n_tokens, ncol = n_rules)
  for (j in seq_len(n_rules)) {
    hits <- stringr::str_detect(match_text, rules$regex[j])
    hits[is.na(hits)] <- FALSE
    match_mat[, j] <- hits
  }

  # Stage 2: per-token best rule per competition group (highest priority wins
  # within each group; groups fire independently)
  priorities <- rules$priority
  groups <- rules$competition_group

  token_best_rules <- lapply(seq_len(n_tokens), function(i) {
    matched <- which(match_mat[i, ])
    if (length(matched) == 0L) {
      return(integer(0L))
    }
    vapply(
      split(matched, groups[matched]),
      function(grp_idx) grp_idx[which.max(priorities[grp_idx])],
      integer(1L)
    )
  })

  has_match <- lengths(token_best_rules) > 0L
  n_matched_rules <- lengths(token_best_rules[has_match])
  rule_indices <- unlist(token_best_rules[has_match], use.names = FALSE)

  token_matches <- tibble::tibble(
    .pk_row_id = rep(tokens_tbl$.pk_row_id[has_match], times = n_matched_rules),
    clone_id = rep(tokens_tbl$clone_id[has_match], times = n_matched_rules),
    aberr_norm = rep(tokens_tbl$aberr_norm[has_match], times = n_matched_rules),
    flag_name = rules$flag_name[rule_indices],
    priority = priorities[rule_indices]
  )

  # Resolve priority per (sample, clone, aberr_norm, flag_name) group
  resolved <- token_matches |>
    dplyr::group_by(.pk_row_id, clone_id, aberr_norm, flag_name) |>
    dplyr::arrange(dplyr::desc(priority)) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  # Stage 3: pivot to sample-level flags
  flags_by_sample <- resolved |>
    dplyr::distinct(.pk_row_id, flag_name) |>
    dplyr::mutate(value = 1L) |>
    tidyr::pivot_wider(
      names_from = flag_name,
      values_from = value,
      values_fill = 0L
    )

  if (nrow(flags_by_sample) == 0) {
    flags_by_sample <- tibble::tibble(
      .pk_row_id = unique(tokens_tbl$.pk_row_id)
    )
  }
  for (nm in rule_flag_names) {
    if (!nm %in% names(flags_by_sample)) flags_by_sample[[nm]] <- 0L
  }
  flags_by_sample
}

compute_aneuploidy <- function(tokens_tbl, chroms) {
  mono_cols <- paste0("mono", chroms)
  tris_cols <- paste0("tris", chroms)

  aneuploidy_tbl <- tokens_tbl |>
    dplyr::filter(aberr_raw != "") |>
    dplyr::transmute(
      .pk_row_id,
      sign = dplyr::case_when(
        stringr::str_detect(aberr_raw, "^-") ~ "-",
        stringr::str_detect(aberr_raw, "^\\+") ~ "+",
        TRUE ~ ""
      ),
      chrom = stringr::str_extract(aberr_raw, "(?<=^[-+])[0-9XY]+"),
      is_mono_tri = sign %in% c("-", "+") & !is.na(chrom),
      prefix = dplyr::case_when(
        sign == "-" ~ "mono",
        sign == "+" ~ "tris",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(is_mono_tri) |>
    dplyr::mutate(flag = paste0(prefix, chrom), value = 1L) |>
    dplyr::distinct(.pk_row_id, flag, .keep_all = TRUE) |>
    dplyr::select(.pk_row_id, flag, value) |>
    tidyr::pivot_wider(names_from = flag, values_from = value, values_fill = 0L)

  for (nm in c(mono_cols, tris_cols)) {
    if (!nm %in% names(aneuploidy_tbl)) aneuploidy_tbl[[nm]] <- 0L
  }
  aneuploidy_tbl
}

compute_comma_counts <- function(tokens_tbl, verbose = FALSE) {
  sex_first_regex <- paste0(
    "^(",
    paste(.sex_complements, collapse = "|"),
    ")$"
  )

  clone_counts <- tokens_tbl |>
    dplyr::filter(!is.na(clone_id)) |>
    dplyr::group_by(.pk_row_id, clone_id) |>
    dplyr::summarise(
      tokens_nonempty = sum(aberr_raw != ""),
      first_token = dplyr::first(aberr_raw[aberr_raw != ""]),
      first_is_sex = !is.na(first_token) &
        stringr::str_detect(first_token, sex_first_regex),
      has_idem = any(aberr_raw == "idem", na.rm = TRUE),
      raw_comma_count = as.integer(
        tokens_nonempty - ifelse(first_is_sex, 1L, 0L)
      ),
      .groups = "drop"
    )

  stemline_counts <- clone_counts |>
    dplyr::filter(clone_id == 1) |>
    dplyr::select(.pk_row_id, stemline_aberr_count = raw_comma_count)

  clone_counts <- clone_counts |>
    dplyr::left_join(stemline_counts, by = ".pk_row_id") |>
    dplyr::mutate(
      idem_invalid = has_idem &
        (clone_id == 1 | is.na(stemline_aberr_count))
    )

  idem_invalid_samples <- clone_counts |>
    dplyr::filter(idem_invalid) |>
    dplyr::distinct(.pk_row_id)

  if (nrow(idem_invalid_samples) > 0 && isTRUE(verbose)) {
    warning(
      sprintf(
        "Found 'idem' without valid stemline in %d sample(s) (row indices: %s). 'idem' in clone 1 or without a preceding clone cannot be expanded.",
        nrow(idem_invalid_samples),
        paste(idem_invalid_samples$.pk_row_id, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  clone_counts |>
    dplyr::mutate(
      comma_count_aberrations = dplyr::case_when(
        has_idem &
          clone_id > 1 &
          !is.na(stemline_aberr_count) ~ raw_comma_count +
          stemline_aberr_count -
          1L,
        TRUE ~ raw_comma_count
      )
    ) |>
    dplyr::group_by(.pk_row_id) |>
    dplyr::summarise(
      comma_count_aberrations = max(tidyr::replace_na(
        comma_count_aberrations,
        0L
      )),
      .groups = "drop"
    )
}

compute_unique_counts <- function(tokens_tbl) {
  tokens_tbl |>
    dplyr::filter(
      !(aberr_norm %in%
        c("", "chromosome_count_range", "idem", "sl", .sex_complements))
    ) |>
    dplyr::group_by(.pk_row_id) |>
    dplyr::summarise(
      n_unique_aberr = dplyr::n_distinct(aberr_norm),
      .groups = "drop"
    )
}

blank_rows <- function(original_karyotypes, all_output_cols) {
  n <- length(original_karyotypes)
  out <- tibble::tibble(original_karyotype = original_karyotypes)
  out$preprocessed_karyotype <- NA_character_
  out$ploidy_category <- NA_character_
  other_cols <- setdiff(
    all_output_cols,
    c("original_karyotype", "preprocessed_karyotype", "ploidy_category")
  )
  for (nm in other_cols) {
    out[[nm]] <- rep(NA_integer_, n)
  }
  out
}

#' Parse ISCN Karyotype Strings
#'
#' Parses ISCN karyotype notation into structured binary features for analysis.
#' Extracts specific translocations, deletions, monosomies, trisomies, and other
#' chromosomal aberrations according to configurable rules.
#'
#' @param karyotypes A character vector of karyotype strings, a data frame
#'   containing a karyotype column, or a `karyo_preprocessed` tibble returned
#'   by `preprocess_karyo()`. When a `karyo_preprocessed` object is passed, the
#'   `preprocessed` column is used automatically and cached issue indices and id
#'   column are reused — no additional arguments required.
#' @param rules A data.frame of parsing rules (default: rules_table()). Must contain
#'   columns: flag_name, regex, category, priority, counts_for_monosomal,
#'   competition_group. Use custom data.frame for specialized parsing needs.
#' @param karyotype_column Character. Name of the karyotype column when input is
#'   a plain data frame. If `NULL` (default), auto-detected from common names
#'   (`karyotype`, `iscn`, etc.); in an interactive session the detected column
#'   is shown with a preview and confirmed before use. Ignored for character
#'   vector or `karyo_preprocessed` input.
#' @param id_column Character. Name of the id column when input is a data frame.
#'   If `NULL` (default), auto-detected from common names (`sample_id`,
#'   `patient_id`, `id`, `mrn`, etc.); in an interactive session the detected
#'   column is confirmed before use, with an option to skip it. The id column is
#'   placed first in the output. For `karyo_preprocessed` input, the id is
#'   propagated automatically
#'   from upstream pipeline steps; pass `id_column` explicitly only to override.
#' @param .return Character string specifying return type: `"tibble"` (default)
#'   or `"data.frame"`.
#' @param verbose Logical. If `TRUE`, prints column detection and parsing
#'   summary messages. Default `FALSE`.
#' @param on_issues Character string specifying how to handle any karyotype
#'   issues detected by `check_karyo()`. Three classes of issues are handled:
#'   - **Dirty markers** (e.g. leading dots, HTML entities, missing sex comma):
#'     formatting artifacts that `preprocess_karyo()` can fix automatically.
#'   - **Chimeric separators** (`//`): karyotypes containing independent cell
#'     populations. Under `"fix"`, only the portion before the first `//` is
#'     kept (the dominant clone). Information about secondary clones is lost.
#'   - **Structural errors** (e.g. no chromosome count, `Updated ISCN` marker):
#'     unfixable — these rows always return NA regardless of `on_issues`.
#'
#'   Values:
#'   - `"fix"` (default): Apply `preprocess_karyo()` to dirty rows; truncate
#'     chimeric rows to the first clone. Rows that still have issues after
#'     these corrections are returned as NA. A message is printed summarising
#'     how many rows were fixed and how many could not be fixed.
#'   - `"warn"`: Return NA for all issue rows and emit a message. No fixing
#'     is attempted.
#'   - `"stop"`: Raise an error immediately if any issues are found. No fixing
#'     is attempted.
#'
#' @return A tibble or data.frame with columns:
#'   - original_karyotype: Input karyotype string (always the raw input)
#'   - preprocessed_karyotype: `normalize_iscn()` output of the karyotype string,
#'     consistent across all `on_issues` modes. In `"fix"` mode this is also
#'     dirty-fixed and chimeric-truncated; in `"warn"`/`"stop"` modes only
#'     `normalize_iscn()` is applied. `NA_character_` for issue rows.
#'   - ploidy_category: Classification (diploid, hyperdiploid, etc.) based on most abnormal clone
#'   - chromosome_count: Integer chromosome count extracted from the karyotype
#'   - One column per aberration flag (0/1 binary)
#'   - monoX, monoY, mono1-22: Monosomy flags for each chromosome
#'   - trisX, trisY, tris1-22: Trisomy flags for each chromosome
#'   - normal_karyotype: 1 if 46,XX or 46,XY, else 0
#'   - total_metaphases: Count from bracket notation
#'   - comma_count_aberrations: Number of comma-separated aberrations
#'   - complex_karyotype: 1 if >=3 unique aberrations, else 0
#'   - monosomal_karyotype: 1 if meets monosomal criteria, else 0
#'   - mixed_ploidy: 1 if clones have different ploidy categories, else 0
#'   - fixable_error: 1 if the row had at least one fixable issue (dirty marker
#'     or chimeric separator) and no unfixable issue, regardless of whether the
#'     fixable issue was auto-corrected. 0 when `unfixable_error = 1`.
#'   - unfixable_error: 1 if row had unfixable structural issues (row is NA), else 0
#'   - chimeric_karyotype: 1 if row contained a `//` chimeric separator (regular
#'     chimeric rows are truncated to the host clone; `zero_host_chimera` rows
#'     are also flagged here and returned as NA via `unfixable_error`)
#'
#' @examples
#' # Basic usage with character vector
#' karyotypes <- c("46,XX", "47,XY,+21", "46,XX,t(15;17)(q24;q21)")
#' result <- parse_karyo(karyotypes)
#'
#' # Usage with data.frame
#' df <- data.frame(
#'   sample_id = c("S1", "S2", "S3"),
#'   iscn = c("46,XX", "47,XY,+21", "46,XX,t(15;17)(q24;q21)")
#' )
#' result <- parse_karyo(df, karyotype_column = "iscn")
#'
#' # Return as data.frame instead of tibble
#' result <- parse_karyo(df, karyotype_column = "iscn", .return = "data.frame")
#'
#' # Use verbose mode for debugging
#' result <- parse_karyo(df, karyotype_column = "iscn", verbose = TRUE)
#'
#' @export
parse_karyo <- function(
  karyotypes,
  rules = rules_table(),
  karyotype_column = NULL,
  id_column = NULL,
  .return = c("tibble", "data.frame"),
  verbose = FALSE,
  on_issues = c("fix", "warn", "stop")
) {
  .return <- match.arg(.return)
  on_issues <- match.arg(on_issues)

  # Validate rules early (needed for empty result structure) ------------------
  required <- c(
    "flag_name",
    "regex",
    "category",
    "priority",
    "counts_for_monosomal",
    "competition_group"
  )
  if (!all(required %in% names(rules))) {
    stop(
      "Rules missing required columns: ",
      paste(setdiff(required, names(rules)), collapse = ", ")
    )
  }
  rule_flag_names <- unique(rules$flag_name)

  # Helper to create empty result with correct structure ----------------------
  chroms <- c(as.character(1:22), "X", "Y")
  all_output_cols <- c(
    "original_karyotype",
    "preprocessed_karyotype",
    "ploidy_category",
    "chromosome_count",
    rule_flag_names,
    paste0("mono", chroms),
    paste0("tris", chroms),
    "normal_karyotype",
    "total_metaphases",
    "comma_count_aberrations",
    "complex_karyotype",
    "monosomal_karyotype",
    "mixed_ploidy",
    "fixable_error",
    "unfixable_error",
    "chimeric_karyotype"
  )

  empty_result <- function() {
    out <- tibble::tibble(original_karyotype = character())
    out$preprocessed_karyotype <- character()
    out$ploidy_category <- character()
    for (nm in setdiff(
      all_output_cols,
      c("original_karyotype", "preprocessed_karyotype", "ploidy_category")
    )) {
      out[[nm]] <- integer()
    }
    attr(out, "karyoparser_version") <- .karyoparser_version
    if (.return == "tibble") out else as.data.frame(out)
  }

  # ---- Input routing: extract raw_vec, original_vec, id info ----------------
  # karyo_preprocessed detection MUST happen before raw_vec is set so the
  # "preprocessed" column is used automatically (it is not in
  # .karyotype_col_candidates).
  is_preprocessed <- is.data.frame(karyotypes) &&
    inherits(karyotypes, "karyo_preprocessed")

  id_values <- NULL
  id_col_name <- NULL

  if (is_preprocessed) {
    raw_vec <- as.character(karyotypes[["preprocessed"]])
    original_vec <- as.character(karyotypes[["original"]])
    # id: explicit id_column param overrides cached attrs
    if (!is.null(id_column)) {
      if (!id_column %in% names(karyotypes)) {
        stop(
          sprintf(
            "ID column '%s' not found in karyo_preprocessed input.",
            id_column
          ),
          call. = FALSE
        )
      }
      id_col_name <- id_column
      id_values <- as.character(karyotypes[[id_column]])
    } else {
      id_col_name <- attr(karyotypes, ".kp_id_col")
      id_values <- attr(karyotypes, ".kp_id_values")
    }
  } else if (is.data.frame(karyotypes)) {
    extracted <- .extract_df_input(
      karyotypes,
      karyotype_column,
      id_column,
      verbose,
      "parse_karyo"
    )
    raw_vec <- extracted$raw_vec
    original_vec <- raw_vec
    id_col_name <- extracted$id_col_name
    id_values <- extracted$id_values
    karyotype_column <- extracted$karyotype_col_name
  } else {
    raw_vec <- as.character(karyotypes)
    original_vec <- raw_vec
  }

  if (length(raw_vec) > 0) {
    n_total_vec <- length(raw_vec)

    if (is_preprocessed) {
      # Reuse cached assessment indices — skip .assess_karyotypes() entirely
      fixable_row_indices <- attr(karyotypes, ".kp_fixable_rows")
      if (is.null(fixable_row_indices)) {
        fixable_row_indices <- integer(0)
      }
      unfixable_row_indices <- attr(karyotypes, ".kp_unfixable_rows")
      if (is.null(unfixable_row_indices)) {
        unfixable_row_indices <- integer(0)
      }
      chimeric_all_indices <- attr(karyotypes, ".kp_chimeric_all_rows")
      if (is.null(chimeric_all_indices)) {
        chimeric_all_indices <- integer(0)
      }
      chimeric_rows_cached <- attr(karyotypes, ".kp_chimeric_rows")
      if (is.null(chimeric_rows_cached)) {
        chimeric_rows_cached <- integer(0)
      }

      if (on_issues == "stop") {
        na_indices <- which(is.na(raw_vec))
        if (length(na_indices) > 0) {
          stop(
            sprintf(
              "%d of %d karyotype(s) are unfixable (NA in preprocessed input). Use on_issues = \"fix\" or \"warn\".",
              length(na_indices),
              n_total_vec
            ),
            call. = FALSE
          )
        }
        issue_row_indices <- integer(0)
        fixable_row_indices <- integer(0)
        unfixable_row_indices <- integer(0)
        chimeric_all_indices <- integer(0)
      } else {
        # "fix" and "warn": NA rows in preprocessed column are already the issue rows
        issue_row_indices <- which(is.na(raw_vec))

        if (isTRUE(verbose)) {
          n_final_issues <- length(unfixable_row_indices)
          n_fixed <- length(setdiff(fixable_row_indices, unfixable_row_indices))
          n_chimeric_parsed <- length(
            setdiff(chimeric_rows_cached, unfixable_row_indices)
          )
          if (n_fixed == 0 && n_final_issues == 0) {
            message(sprintf("Parsed %d karyotype(s).", n_total_vec))
          } else {
            parts <- character(0)
            if (n_fixed > 0) {
              parts <- c(parts, sprintf("%d pre-fixed", n_fixed))
            }
            if (n_final_issues > 0) {
              parts <- c(parts, sprintf("%d as NA", n_final_issues))
            }
            message(sprintf(
              "Parsed %d/%d karyotype(s). %s.",
              n_total_vec - n_final_issues,
              n_total_vec,
              paste(parts, collapse = ", ")
            ))
            if (n_chimeric_parsed > 0) {
              message(sprintf(
                "  Chimeric: %d \u2014 host clone extracted, donor discarded.",
                n_chimeric_parsed
              ))
            }
          }
        }
      }
    } else {
      # Standard path: run unified assessment
      assessment <- .assess_karyotypes(raw_vec)

      fixable_row_indices <- unique(
        union(assessment$dirty_row_indices, assessment$chimeric_row_indices)
      )
      unfixable_row_indices <- assessment$unfixable_row_indices
      zhc_indices <- unique(assessment$reported_issues$row_index[
        assessment$reported_issues$issue_type == "zero_host_chimera"
      ])
      chimeric_all_indices <- unique(c(
        assessment$chimeric_row_indices,
        zhc_indices
      ))

      if (on_issues == "stop") {
        if (nrow(assessment$reported_issues) > 0) {
          issue_counts <- assessment$reported_issues |>
            dplyr::count(issue_type, name = "n") |>
            dplyr::arrange(dplyr::desc(n))
          stop(
            sprintf(
              "%d of %d karyotype(s) have issues (%s). Use on_issues='fix' or 'warn', or call check_karyo() for details.",
              length(unique(assessment$reported_issues$row_index)),
              n_total_vec,
              paste(
                paste0(issue_counts$issue_type, ": ", issue_counts$n),
                collapse = ", "
              )
            ),
            call. = FALSE
          )
        }
        raw_vec <- normalize_iscn(raw_vec)
        issue_row_indices <- integer(0)
        fixable_row_indices <- integer(0)
        unfixable_row_indices <- integer(0)
        chimeric_all_indices <- integer(0)
      } else if (on_issues == "fix") {
        raw_vec <- assessment$processed
        issue_row_indices <- unfixable_row_indices

        if (isTRUE(verbose)) {
          n_chimeric_parsed <- length(setdiff(
            assessment$chimeric_row_indices,
            assessment$unfixable_row_indices
          ))
          n_fixed <- length(setdiff(
            fixable_row_indices,
            assessment$unfixable_row_indices
          ))
          n_final_issues <- length(unfixable_row_indices)
          if (n_fixed == 0 && n_final_issues == 0) {
            message(sprintf("Parsed %d karyotype(s).", n_total_vec))
          } else {
            parts <- character(0)
            if (n_fixed > 0) {
              parts <- c(parts, sprintf("%d fixed", n_fixed))
            }
            if (n_final_issues > 0) {
              parts <- c(parts, sprintf("%d as NA", n_final_issues))
            }
            message(sprintf(
              "Parsed %d/%d karyotype(s). %s.",
              n_total_vec - n_final_issues,
              n_total_vec,
              paste(parts, collapse = ", ")
            ))
            if (n_chimeric_parsed > 0) {
              message(sprintf(
                "  Chimeric: %d \u2014 host clone extracted, donor discarded.",
                n_chimeric_parsed
              ))
            }
          }
        }
      } else {
        # "warn": detect issues on raw, no fixing; all issue rows become NA
        raw_vec <- normalize_iscn(raw_vec)
        issue_row_indices <- unique(assessment$reported_issues$row_index)

        if (isTRUE(verbose)) {
          n_final_issues <- length(issue_row_indices)
          if (n_final_issues == 0) {
            message(sprintf("Parsed %d karyotype(s).", n_total_vec))
          } else {
            message(sprintf(
              "%d of %d karyotype(s) had issues and were returned as NA.",
              n_final_issues,
              n_total_vec
            ))
          }
        }
      }
    }
  } else {
    issue_row_indices <- integer(0)
    fixable_row_indices <- integer(0)
    unfixable_row_indices <- integer(0)
    chimeric_all_indices <- integer(0)
  }

  # Ingest --------------------------------------------------------------------
  # raw_vec, original_vec, id_values, id_col_name were all resolved in the
  # input-routing block above. Just build input_df here.
  input_df <- tibble::tibble(
    .pk_row_id = seq_along(raw_vec),
    original_karyotype = original_vec,
    preprocessed_karyotype = raw_vec
  )

  # Apply issue-row filtering -------------------------------------------------
  if (length(issue_row_indices) > 0) {
    issue_rows_df <- input_df |>
      dplyr::filter(.pk_row_id %in% issue_row_indices)
    input_df <- input_df |>
      dplyr::filter(!(.pk_row_id %in% issue_row_indices))
  } else {
    issue_rows_df <- input_df[0, ]
  }

  # Check for empty input -----------------------------------------------------
  if (nrow(input_df) == 0) {
    if (nrow(issue_rows_df) > 0) {
      # All rows had issues - return blank rows with issues attached
      out <- blank_rows(issue_rows_df$original_karyotype, all_output_cols)
      out$.pk_row_id <- issue_rows_df$.pk_row_id
      out$fixable_error <- as.integer(
        out$.pk_row_id %in%
          fixable_row_indices &
          !out$.pk_row_id %in% unfixable_row_indices
      )
      out$unfixable_error <- as.integer(
        out$.pk_row_id %in% unfixable_row_indices
      )
      out$chimeric_karyotype <- as.integer(
        out$.pk_row_id %in% chimeric_all_indices
      )
      if (!is.null(id_values)) {
        out[[id_col_name]] <- id_values[out$.pk_row_id]
        out <- out |> dplyr::relocate(dplyr::all_of(id_col_name), .before = 1)
      }
      out <- out |>
        dplyr::select(-.pk_row_id) |>
        dplyr::relocate(preprocessed_karyotype, .after = original_karyotype)
      attr(out, "karyoparser_version") <- .karyoparser_version
      return(
        if (.return == "tibble") tibble::as_tibble(out) else as.data.frame(out)
      )
    }
    if (isTRUE(verbose)) {
      message("Input is empty. Returning empty result.")
    }
    return(empty_result())
  }

  # ---- Deduplication ---------------------------------------------------------
  n_total <- nrow(input_df)
  unique_karyotypes <- unique(input_df$preprocessed_karyotype)
  n_unique <- length(unique_karyotypes)
  deduped <- n_unique < n_total

  if (deduped) {
    if (isTRUE(verbose)) {
      message(sprintf(
        "Parsing %d unique karyotypes (%d total rows).",
        n_unique,
        n_total
      ))
    }
    # Parse unique normalized strings only, then join back
    dedup_df <- tibble::tibble(
      .pk_row_id = seq_along(unique_karyotypes),
      original_karyotype = unique_karyotypes
    )
  } else {
    if (isTRUE(verbose)) {
      message(sprintf("Parsing %d karyotype(s).", n_total))
    }
    dedup_df <- tibble::tibble(
      .pk_row_id = input_df$.pk_row_id,
      original_karyotype = input_df$preprocessed_karyotype
    )
  }

  # ---- Parsing pipeline -------------------------------------------------------
  sample_meta <- build_sample_meta(dedup_df, verbose)
  tokens <- build_clone_tokens(sample_meta)
  flags <- match_rules(tokens, rules, rule_flag_names)
  aneuploidy <- compute_aneuploidy(tokens, chroms)
  comma_counts <- compute_comma_counts(tokens, verbose)
  unique_aberr <- compute_unique_counts(tokens)

  # ---- Assembly --------------------------------------------------------------
  mono_cols <- paste0("mono", chroms)
  tris_cols <- paste0("tris", chroms)

  all_flags <- dedup_df |>
    dplyr::select(.pk_row_id) |>
    dplyr::left_join(flags, by = ".pk_row_id") |>
    dplyr::left_join(aneuploidy, by = ".pk_row_id") |>
    dplyr::mutate(dplyr::across(-.pk_row_id, ~ tidyr::replace_na(., 0L)))

  structural_flags <- rules |>
    dplyr::filter(counts_for_monosomal) |>
    dplyr::pull(flag_name) |>
    unique()

  autosomal_mono_cols <- paste0("mono", as.character(1:22))
  available_mono_cols <- intersect(autosomal_mono_cols, names(all_flags))
  available_struct_flags <- intersect(structural_flags, names(all_flags))

  abnormality_names <- c(
    rule_flag_names,
    mono_cols,
    tris_cols,
    "normal_karyotype",
    "comma_count_aberrations",
    "complex_karyotype",
    "monosomal_karyotype",
    "mixed_ploidy"
  )

  result <- all_flags |>
    dplyr::mutate(
      autosomal_monosomies_sample = if (length(available_mono_cols) > 0) {
        rowSums(
          dplyr::across(dplyr::all_of(available_mono_cols), as.integer),
          na.rm = TRUE
        )
      } else {
        0L
      },
      structural_aberrations_sample = if (length(available_struct_flags) > 0) {
        rowSums(
          dplyr::across(dplyr::all_of(available_struct_flags), as.integer),
          na.rm = TRUE
        )
      } else {
        0L
      },
      monosomal_karyotype = as.integer(
        autosomal_monosomies_sample >= 2 |
          (autosomal_monosomies_sample >= 1 &
            structural_aberrations_sample >= 1)
      )
    ) |>
    dplyr::select(
      -autosomal_monosomies_sample,
      -structural_aberrations_sample
    )

  # CBF-AML override: cases with t(8;21), inv(16)(p13q22), or t(16;16) are
  # never monosomal per clinical guidelines, regardless of co-firing rules.
  cbf_flags <- c("t(8;21)(q22;q22)", "inv(16)(p13q22)", "t(16;16)(p13;q22)")
  available_cbf <- intersect(cbf_flags, names(result))
  if (length(available_cbf) > 0) {
    result <- result |>
      dplyr::mutate(
        monosomal_karyotype = dplyr::if_else(
          rowSums(
            dplyr::across(dplyr::all_of(available_cbf), as.integer),
            na.rm = TRUE
          ) >=
            1L,
          0L,
          monosomal_karyotype
        )
      )
  }

  result <- result |>
    dplyr::left_join(
      unique_aberr |>
        dplyr::transmute(
          .pk_row_id,
          complex_karyotype = as.integer(n_unique_aberr >= 3)
        ),
      by = ".pk_row_id"
    ) |>
    dplyr::left_join(comma_counts, by = ".pk_row_id") |>
    dplyr::left_join(
      sample_meta |>
        dplyr::select(
          .pk_row_id,
          normal_karyotype,
          ploidy_category,
          chromosome_count,
          total_metaphases,
          mixed_ploidy
        ),
      by = ".pk_row_id"
    )

  for (nm in abnormality_names) {
    if (!nm %in% names(result)) result[[nm]] <- 0L
  }

  cols_to_convert <- intersect(abnormality_names, names(result))
  result <- result |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(cols_to_convert),
        ~ as.integer(tidyr::replace_na(., 0L))
      ),
      ploidy_category = tidyr::replace_na(ploidy_category, "unknown")
    ) |>
    dplyr::left_join(dedup_df, by = ".pk_row_id") |>
    dplyr::rename(preprocessed_karyotype = original_karyotype) |>
    dplyr::relocate(preprocessed_karyotype, .before = 1) |>
    dplyr::relocate(ploidy_category, .after = preprocessed_karyotype) |>
    dplyr::relocate(chromosome_count, .after = ploidy_category)

  keep_cols <- intersect(
    c(
      "preprocessed_karyotype",
      "ploidy_category",
      "chromosome_count",
      abnormality_names,
      "total_metaphases"
    ),
    names(result)
  )

  # ---- Join back to original rows if deduped ---------------------------------
  if (deduped) {
    parsed_unique <- result[, keep_cols]
    valid_out <- input_df |>
      dplyr::left_join(parsed_unique, by = "preprocessed_karyotype")
    valid_out <- valid_out[, c(".pk_row_id", "original_karyotype", keep_cols)]
  } else {
    valid_out <- dplyr::bind_cols(
      tibble::tibble(
        .pk_row_id = input_df$.pk_row_id,
        original_karyotype = input_df$original_karyotype
      ),
      result[, keep_cols]
    )
  }

  # Recombine with issue rows --------------------------------------------------
  if (nrow(issue_rows_df) > 0) {
    br <- blank_rows(issue_rows_df$original_karyotype, all_output_cols)
    br$.pk_row_id <- issue_rows_df$.pk_row_id
    out <- dplyr::bind_rows(valid_out, br) |>
      dplyr::arrange(.pk_row_id)
  } else {
    out <- valid_out
  }

  # Error flag columns ---------------------------------------------------------
  out$fixable_error <- as.integer(
    out$.pk_row_id %in%
      fixable_row_indices &
      !out$.pk_row_id %in% unfixable_row_indices
  )
  out$unfixable_error <- as.integer(out$.pk_row_id %in% unfixable_row_indices)
  out$chimeric_karyotype <- as.integer(out$.pk_row_id %in% chimeric_all_indices)

  # Attach ID column if available ---------------------------------------------
  if (!is.null(id_values)) {
    out[[id_col_name]] <- id_values[out$.pk_row_id]
    out <- out |> dplyr::relocate(dplyr::all_of(id_col_name), .before = 1)
  }

  out <- out |>
    dplyr::select(-.pk_row_id) |>
    dplyr::relocate(preprocessed_karyotype, .after = original_karyotype)

  # Provenance
  attr(out, "karyoparser_version") <- .karyoparser_version
  if (.return == "tibble") tibble::as_tibble(out) else as.data.frame(out)
}
