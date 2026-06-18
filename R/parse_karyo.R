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
#'   column are reused -- no additional arguments required.
#' @param rules A `karyo_rules` object (default: [myeloid_rules]). Use
#'   [validate_rules()] to validate and convert a custom data frame into an
#'   accepted rules object.
#' @param karyotype_column Character. Name of the karyotype column when input is
#'   a plain data frame. If `NULL` (default), auto-detected from common names
#'   (`karyotype`, `iscn`, etc.). Ignored for character vector or
#'   `karyo_preprocessed` input.
#' @param id_column Character. Name of the id column when input is a data frame.
#'   If `NULL` (default), auto-detected from common names (`sample_id`,
#'   `patient_id`, `id`, `mrn`, etc.) and used silently. The id column is
#'   placed first in the output. For `karyo_preprocessed` input, the id is
#'   propagated automatically
#'   from upstream pipeline steps; pass `id_column` explicitly only to override.
#' @param verbose Logical. If `TRUE`, prints column detection and parsing
#'   summary messages. Default `FALSE`.
#' @param on_issues Character string specifying how to handle **dirty markers**
#'   (formatting artifacts such as leading dots, HTML entities, missing sex
#'   comma) and **structural errors** (e.g. no chromosome count, `Updated ISCN`
#'   marker). Chimeric handling is independent of this argument -- see
#'   `on_chimeric`. Structural errors are always unfixable and return NA
#'   regardless of `on_issues`.
#'
#'   Values:
#'   - `"stop"` (default): Raise an error immediately if any non-chimeric issue
#'     is found, listing fixable and unfixable issue types and counts, and
#'     suggesting next steps. No fixing is attempted. Chimeric rows never
#'     trigger the error -- they are clone-selected and parsed per
#'     `on_chimeric`. **Exception**: when input is a `karyo_preprocessed`
#'     object, `"stop"` does not error -- unfixable rows are returned as NA and
#'     a warning is emitted.
#'   - `"preprocess"`: Apply `preprocess_karyo()` to dirty rows. Rows that still
#'     have issues after these corrections are returned as NA. A message is
#'     printed summarising how many rows were fixed and how many could not be
#'     fixed.
#'   - `"warn"`: Return NA for all dirty/structural issue rows and emit a
#'     warning. No fixing is attempted.
#'
#'   Across all three modes, chimeric rows are clone-selected per `on_chimeric`
#'   and parsed; a chimeric row becomes NA only when its selected clone is
#'   itself unparseable (or it has two or more `//` separators).
#' @param on_chimeric Character string controlling which clone is parsed for
#'   chimeric karyotypes (those containing a `//` separator). One of:
#'   - `"default"`: host clone (before `//`) for normal chimeras; donor clone
#'     (after `//`) for zero-host chimeras (`.//` prefix).
#'   - `"host"`: host clone only; zero-host chimeras return NA.
#'   - `"donor"`: everything after `//` (the donor population).
#'   Rows with two or more `//` separators (`multiple_chimeric_separator`) are
#'   always unfixable and return NA. For `karyo_preprocessed` input the clone
#'   selection was baked in at preprocess time and is inherited; passing an
#'   explicit, conflicting `on_chimeric` raises an error.
#'
#' @return A tibble with columns:
#'   - original_karyotype: Input karyotype string (always the raw input)
#'   - preprocessed_karyotype: `normalize_iscn()` output of the
#'     karyotype string, consistent across all `on_issues` modes. In
#'     `"preprocess"` mode this is also dirty-fixed; in `"warn"`/`"stop"` modes
#'     only `normalize_iscn()` is applied (chimeric rows are clone-selected in
#'     every mode). `NA_character_` for issue rows.
#'   - ploidy_category: Classification (diploid, hyperdiploid, etc.)
#'     based on most abnormal clone
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
#'   - balanced_translocation: 1 if the row carries a balanced translocation --
#'     a bare `t(...)` token, or a reciprocal der pair where both partner
#'     chromosomes appear as centromere donors (e.g.
#'     `der(5)t(5;17)...,der(17)t(5;17)...`). Independent of
#'     `derivative_chromosome`. A row may be both balanced and unbalanced.
#'   - unbalanced_translocation: 1 if the row carries an unbalanced
#'     translocation -- a lone `der(a)t(a;b)` (reciprocal der absent) or a
#'     whole-arm `der(a;b)`. The specific fusion flag still fires (e.g.
#'     `der(9)t(9;22)` keeps `t(9;22)(q34;q11) = 1`).
#'   - unbal_loss_<arm>: one column per chromosome arm (`unbal_loss_1p`,
#'     `unbal_loss_1q`, ..., `unbal_loss_22q`, `unbal_loss_Xp`, `unbal_loss_Xq`,
#'     `unbal_loss_Yp`, `unbal_loss_Yq`), set to 1 when an unbalanced der
#'     translocation implies loss of that arm. Kept SEPARATE from
#'     `del(...)`/`mono*` -- an unbalanced-derived 5q loss does not set
#'     `del(5q)`. Derivation needs explicit breakpoints and is limited to simple
#'     single-junction `der(a)t(a;b)`; multi-junction chains and three-way
#'     `t(a;b;c)` derivatives are flagged unbalanced but derive no loss.
#'     Copy-number-aware gains are out of scope.
#'   - unbal_partial_loss: 1 if an unbalanced der translocation implies any
#'     partial loss (OR across all `unbal_loss_*` columns).
#'   - fixable_error: 1 if the row had at least one fixable issue (dirty marker
#'     or chimeric separator) and no unfixable issue, regardless of whether the
#'     fixable issue was auto-corrected. 0 when `unfixable_error = 1`.
#'   - unfixable_error: 1 if row had unfixable structural issues
#'     (row is NA), else 0
#'   - chimeric_karyotype: 1 if row contained a `//` chimeric separator
#'     (including zero-host `.//` chimeras and `multiple_chimeric_separator`
#'     rows), else 0
#'   - chimeric_clone: which clone was parsed for a chimeric row -- `"host"`,
#'     `"donor"`, or `NA_character_` (non-chimeric rows, or chimeric rows with
#'     no selectable clone such as `multiple_chimeric_separator` or a zero-host
#'     chimera under `on_chimeric = "host"`)
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
#' # Use verbose mode for debugging
#' result <- parse_karyo(df, karyotype_column = "iscn", verbose = TRUE)
#'
#' @export
parse_karyo <- function(
  karyotypes,
  rules = myeloid_rules,
  karyotype_column = NULL,
  id_column = NULL,
  verbose = FALSE,
  on_issues = c("stop", "preprocess", "warn"),
  on_chimeric = c("default", "host", "donor")
) {
  on_issues <- match.arg(on_issues)
  on_chimeric_provided <- !missing(on_chimeric)
  on_chimeric <- match.arg(on_chimeric)

  # Validate rules early (needed for empty result structure) ------------------
  if (!inherits(rules, "karyo_rules")) {
    stop(
      "`rules` must be a `karyo_rules` object. ",
      "Pass `myeloid_rules` or use `validate_rules()` to validate a custom rules table.",
      call. = FALSE
    )
  }
  rule_flag_names <- unique(rules$flag_name)
  chroms <- c(as.character(1:22), "X", "Y")

  # Single source of truth for the output columns: a provenance-tagged catalog
  # of every column (see .column_catalog()). `all_output_cols` is the full
  # ordered column list; `abnormality_names` is the subset of binary 0/1 flags
  # that share the parsed-row treatment (backfill missing with 0L,
  # integer-coerce, NA -> 0); `char_cols` is the set of character-typed columns.
  # Adding a feature in .column_catalog() flows to all three automatically.
  catalog <- .column_catalog(rules)
  all_output_cols <- catalog$column
  abnormality_names <- catalog$column[
    catalog$class %in%
      c("rule", "general", "aneuploidy", "summary", "balance", "loss")
  ]
  char_cols <- catalog$column[catalog$type == "character"]

  empty_result <- function() {
    out <- tibble::tibble(.rows = 0L)
    # Build in catalog order so the empty result matches a normal parsed result.
    for (nm in all_output_cols) {
      out[[nm]] <- if (nm %in% char_cols) character() else integer()
    }
    attr(out, "karyoparser_version") <- .karyoparser_version
    out
  }

  # ---- Input routing: extract raw_vec, original_vec, id info ----------------
  # Detection before raw_vec is set: "preprocessed" column is not in .karyotype_col_candidates.
  is_preprocessed <- is.data.frame(karyotypes) &&
    inherits(karyotypes, "karyo_preprocessed")

  id_values <- NULL
  id_col_name <- NULL

  if (is_preprocessed) {
    raw_vec <- as.character(karyotypes[["preprocessed"]])
    original_vec <- as.character(karyotypes[["original"]])
    # Clone selection was baked in at preprocess time; inherit it. An explicit,
    # conflicting on_chimeric cannot be honored without re-deriving from the raw
    # input, so stop and tell the user how to proceed.
    cached_on_chimeric <- attr(karyotypes, ".kp_on_chimeric")
    if (is.null(cached_on_chimeric)) {
      cached_on_chimeric <- "default"
    }
    if (on_chimeric_provided && on_chimeric != cached_on_chimeric) {
      stop(
        sprintf(
          paste(
            "on_chimeric = \"%s\" conflicts with the karyo_preprocessed input,",
            "which was created with on_chimeric = \"%s\".",
            "\nThe donor/host clone selection is baked into the `preprocessed`",
            "column, so it cannot be changed here.",
            "\nRe-run preprocess_karyo(..., on_chimeric = \"%s\"), or pass the",
            "raw karyotypes to parse_karyo() with on_chimeric = \"%s\"."
          ),
          on_chimeric,
          cached_on_chimeric,
          on_chimeric,
          on_chimeric
        ),
        call. = FALSE
      )
    }
    on_chimeric <- cached_on_chimeric
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
      # Reuse cached assessment indices -- skip .assess_karyotypes() entirely.
      # Clone selection was already applied to the `preprocessed` column.
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
      chimeric_clone_vec <- attr(karyotypes, ".kp_chimeric_clone")
      if (is.null(chimeric_clone_vec)) {
        chimeric_clone_vec <- rep(NA_character_, n_total_vec)
      }

      # Override "stop" for karyo_preprocessed: user has already inspected; NA rows warn, don't error.
      issue_row_indices <- which(is.na(raw_vec))
      if (on_issues == "stop" && length(issue_row_indices) > 0) {
        warning(
          sprintf(
            "%d of %d preprocessed karyotype(s) were unfixable and returned as NA.",
            length(issue_row_indices),
            n_total_vec
          ),
          call. = FALSE
        )
      }

      n_final_issues <- length(unfixable_row_indices)
      n_fixed <- length(setdiff(fixable_row_indices, unfixable_row_indices))
      if (isTRUE(verbose)) {
        if (n_fixed > 0) {
          message("  Fixed:      ", n_fixed)
        }
        # Skip verbose message in "stop" mode -- the unconditional warning already covers it.
        if (n_final_issues > 0 && on_issues != "stop") {
          message("  Unfixable:  ", n_final_issues, "  (returned as NA)")
        }
      }
    } else {
      assessment <- .assess_karyotypes(raw_vec, on_chimeric)

      fixable_row_indices <- unique(c(
        assessment$dirty_row_indices,
        assessment$chimeric_row_indices,
        assessment$zhc_row_indices
      ))
      unfixable_row_indices <- assessment$unfixable_row_indices
      chimeric_all_indices <- unique(c(
        assessment$chimeric_row_indices,
        assessment$zhc_row_indices,
        assessment$multi_row_indices
      ))
      chimeric_clone_vec <- assessment$chimeric_clone
      proc <- assessment$processed

      # Chimeric handling is orthogonal to on_issues: chimeric rows are always
      # clone-selected and parsed (NA only when the selected clone is itself
      # unparseable). on_issues governs dirty + structural rows only.
      if (on_issues == "stop") {
        non_chimeric_issues <- assessment$reported_issues |>
          dplyr::filter(!row_index %in% chimeric_all_indices)
        if (nrow(non_chimeric_issues) > 0) {
          issue_counts <- non_chimeric_issues |>
            dplyr::count(issue_type, name = "n") |>
            dplyr::arrange(dplyr::desc(n))
          n_issue_rows <- length(unique(non_chimeric_issues$row_index))
          fixable_counts <- issue_counts[
            issue_counts$issue_type %in% .fixable_issue_types,
            ,
            drop = FALSE
          ]
          unfixable_counts <- issue_counts[
            !issue_counts$issue_type %in% .fixable_issue_types,
            ,
            drop = FALSE
          ]
          msg_lines <- sprintf(
            "parse_karyo() stopped: %d of %d karyotype(s) have data quality issues.",
            n_issue_rows,
            n_total_vec
          )
          if (nrow(fixable_counts) > 0) {
            msg_lines <- c(
              msg_lines,
              paste0(
                "  Fixable:   ",
                paste(
                  paste0(
                    fixable_counts$issue_type,
                    " (",
                    fixable_counts$n,
                    ")"
                  ),
                  collapse = ", "
                )
              )
            )
          }
          if (nrow(unfixable_counts) > 0) {
            msg_lines <- c(
              msg_lines,
              paste0(
                "  Unfixable: ",
                paste(
                  paste0(
                    unfixable_counts$issue_type,
                    " (",
                    unfixable_counts$n,
                    ")"
                  ),
                  collapse = ", "
                )
              )
            )
          }
          msg_lines <- c(
            msg_lines,
            "Rerun with on_issues = \"preprocess\" to auto-correct fixable rows (unfixable rows will be NA).",
            "Rerun with on_issues = \"warn\" to return NA for all issue rows without fixing.",
            "Call check_karyo() for a full per-row quality report."
          )
          stop(paste(msg_lines, collapse = "\n"), call. = FALSE)
        }
        # Non-chimeric rows are clean; chimeric rows use the selected clone.
        raw_vec <- normalize_iscn(raw_vec)
        raw_vec[chimeric_all_indices] <- proc[chimeric_all_indices]
        issue_row_indices <- which(is.na(raw_vec))
      } else if (on_issues == "preprocess") {
        raw_vec <- proc
        issue_row_indices <- unfixable_row_indices

        n_fixed <- length(setdiff(
          fixable_row_indices,
          assessment$unfixable_row_indices
        ))
        n_final_issues <- length(unfixable_row_indices)
        if (isTRUE(verbose)) {
          if (n_fixed > 0) {
            message("  Fixed:      ", n_fixed)
          }
          if (n_final_issues > 0) {
            message("  Unfixable:  ", n_final_issues, "  (returned as NA)")
          }
        }
      } else {
        # warn: dirty + structural rows are returned NA (no fix attempted), but
        # chimeric rows are still clone-selected and parsed.
        norm_vec <- normalize_iscn(raw_vec)
        norm_vec[chimeric_all_indices] <- proc[chimeric_all_indices]
        raw_vec <- norm_vec

        non_chimeric_issue_rows <- setdiff(
          unique(assessment$reported_issues$row_index),
          chimeric_all_indices
        )
        chimeric_na_rows <- intersect(which(is.na(proc)), chimeric_all_indices)
        issue_row_indices <- unique(c(
          non_chimeric_issue_rows,
          chimeric_na_rows
        ))

        n_final_issues <- length(issue_row_indices)
        if (n_final_issues > 0) {
          warning(
            sprintf(
              "%d of %d karyotype(s) had issues and were returned as NA.",
              n_final_issues,
              n_total_vec
            ),
            call. = FALSE
          )
        }
      }
    }

    # Single informational message for chimeric rows (replaces the old warning).
    n_chimeric <- length(chimeric_all_indices)
    if (n_chimeric > 0) {
      message(sprintf(
        "Chimeric: %d (on_chimeric = \"%s\").",
        n_chimeric,
        on_chimeric
      ))
    }
  } else {
    issue_row_indices <- integer(0)
    fixable_row_indices <- integer(0)
    unfixable_row_indices <- integer(0)
    chimeric_all_indices <- integer(0)
    chimeric_clone_vec <- character(0)
  }

  # Ingest --------------------------------------------------------------------
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
      out <- blank_rows(
        issue_rows_df$original_karyotype,
        all_output_cols,
        char_cols
      )
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
      out$chimeric_clone <- chimeric_clone_vec[out$.pk_row_id]
      if (!is.null(id_values)) {
        out[[id_col_name]] <- id_values[out$.pk_row_id]
        out <- out |> dplyr::relocate(dplyr::all_of(id_col_name), .before = 1)
      }
      out <- out |>
        dplyr::select(-.pk_row_id) |>
        dplyr::relocate(preprocessed_karyotype, .after = original_karyotype)
      attr(out, "karyoparser_version") <- .karyoparser_version
      return(tibble::as_tibble(out))
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
  sample_meta <- build_sample_meta(dedup_df)
  tokens <- build_clone_tokens(sample_meta)
  flags <- match_rules(tokens, rules, rule_flag_names)
  aneuploidy <- compute_aneuploidy(tokens, chroms)
  general_flags <- compute_general_flags(tokens)
  comma_counts <- compute_comma_counts(tokens)
  unique_aberr <- compute_unique_counts(tokens)
  der_translocations <- .der_translocations(tokens)
  balance_flags <- classify_translocation_balance(tokens, der_translocations)
  unbal_loss <- derive_unbalanced_loss(der_translocations)

  # ---- Assembly --------------------------------------------------------------
  all_flags <- dedup_df |>
    dplyr::select(.pk_row_id) |>
    dplyr::left_join(flags, by = ".pk_row_id") |>
    dplyr::left_join(aneuploidy, by = ".pk_row_id") |>
    dplyr::left_join(general_flags, by = ".pk_row_id") |>
    dplyr::mutate(dplyr::across(-.pk_row_id, ~ tidyr::replace_na(., 0L)))

  # A row has a structural aberration (for monosomal-karyotype) if any general
  # structural flag fires. The general_* detections cover every structural
  # category universally (translocation, deletion, inversion, addition,
  # dicentric, isodicentric, isochromosome, ring, insertion, duplication,
  # triplication, derivative), so "is there a structural aberration?" no longer
  # depends on per-rule annotation. The lone clinical exception, CBF-AML, is
  # applied as an explicit override below.
  structural_flags <- .general_flags_for_monosomal

  autosomal_mono_cols <- paste0("mono", as.character(1:22))
  available_mono_cols <- intersect(autosomal_mono_cols, names(all_flags))
  available_struct_flags <- intersect(structural_flags, names(all_flags))

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

  # CBF-AML cases are never monosomal per clinical guidelines, regardless of co-firing rules.
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
    dplyr::left_join(balance_flags, by = ".pk_row_id") |>
    dplyr::left_join(unbal_loss, by = ".pk_row_id") |>
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
    br <- blank_rows(
      issue_rows_df$original_karyotype,
      all_output_cols,
      char_cols
    )
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
  out$chimeric_clone <- chimeric_clone_vec[out$.pk_row_id]

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
  tibble::as_tibble(out)
}

# ---- Internal pipeline helpers ------------------------------------------------

build_sample_meta <- function(input_df) {
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

# Each rule fires independently: a flag is 1 for a row if any of the row's
# tokens matches the rule's regex. Rules do not compete -- a single token can
# light up several flags at once (e.g. a specific t(9;11)(p21;q23) plus the
# family flag t(v;11q23)). Mutual exclusivity, where wanted (e.g. a "*_other"
# variant must not fire for its own canonical breakpoints), is expressed in the
# regex itself via negative lookahead, not via a priority ranking.
match_rules <- function(tokens_tbl, rules, rule_flag_names) {
  match_text <- strip_bands(tokens_tbl$aberr_raw)
  ids <- tokens_tbl$.pk_row_id
  out <- tibble::tibble(.pk_row_id = unique(ids))
  for (nm in rule_flag_names) {
    out[[nm]] <- 0L
  }
  for (j in seq_len(nrow(rules))) {
    hit <- stringr::str_detect(match_text, rules$regex[j])
    hit[is.na(hit)] <- FALSE
    nm <- rules$flag_name[j]
    out[[nm]] <- pmax(
      out[[nm]],
      as.integer(out$.pk_row_id %in% unique(ids[hit]))
    )
  }
  out
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

# General, disease-agnostic structural-aberration detections. Unlike the regex
# rules these are universal (always computed regardless of the rule set) and
# fire INDEPENDENTLY of any specific rule matching the same token -- e.g. a
# t(8;21) token sets both the specific t(8;21)(q22;q22) flag and
# general_translocation. Patterns mirror the historical "general" rules and are
# matched against the same target as match_rules (strip_bands(aberr_raw)).
# Surfaced as the `general_*` output columns (catalog class "general"); kept as
# a named vector so the column set and the catalog stay in sync. To later make
# their surfacing toggleable, gate this group via the catalog class.
.general_flag_patterns <- c(
  general_dicentric = "dic\\(",
  general_isodicentric = "idic\\(",
  general_isochromosome = "i\\(",
  general_pseudodicentric = "psu dic\\(",
  general_ring = "\\br\\(",
  general_insertion = "ins\\(",
  general_duplication = "dup\\(",
  general_triplication = "trp\\(",
  general_translocation = "^[?~]?t\\([0-9XY]",
  general_addition = "add\\(",
  general_inversion = "inv\\(",
  general_deletion = "del\\(",
  general_marker = "\\bmar\\b",
  general_derivative = "^\\+?i?der\\("
)

# General structural flags that count as a structural aberration for the
# monosomal-karyotype rule (everything except a lone marker).
.general_flags_for_monosomal <- setdiff(
  names(.general_flag_patterns),
  "general_marker"
)

compute_general_flags <- function(tokens_tbl) {
  match_text <- strip_bands(tokens_tbl$aberr_raw)
  ids <- tokens_tbl$.pk_row_id
  out <- tibble::tibble(.pk_row_id = unique(ids))
  for (nm in names(.general_flag_patterns)) {
    hit <- stringr::str_detect(match_text, .general_flag_patterns[[nm]])
    hit[is.na(hit)] <- FALSE
    out[[nm]] <- as.integer(out$.pk_row_id %in% unique(ids[hit]))
  }
  out
}

compute_comma_counts <- function(tokens_tbl) {
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
    dplyr::left_join(stemline_counts, by = ".pk_row_id")

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

# All chromosome arms (autosomes 1-22 plus X/Y, each p and q) for which
# unbalanced-derived partial losses get their own `unbal_loss_<arm>` column.
.unbal_loss_arms <- paste0(
  rep(c(as.character(1:22), "X", "Y"), each = 2),
  c("p", "q")
)

# Extract der(a)t(a;b) tokens with a canonical translocation signature and a
# per-(row, clone, signature) `balanced_pair` flag (TRUE when both partner
# chromosomes appear as centromere donors -- i.e. a reciprocal der pair).
# Shared by classify_translocation_balance() and derive_unbalanced_loss().
.der_translocations <- function(tokens_tbl) {
  der_t <- tokens_tbl |>
    dplyr::filter(!is.na(clone_id), aberr_raw != "") |>
    dplyr::mutate(
      raw = strip_bands(aberr_raw),
      is_der = stringr::str_detect(raw, "^\\+?i?der\\("),
      der_inner = stringr::str_match(raw, "^\\+?i?der\\(([^)]+)\\)")[, 2],
      t_c1 = stringr::str_match(raw, "t\\(([0-9XY]+);([0-9XY]+)\\)")[, 2],
      t_c2 = stringr::str_match(raw, "t\\(([0-9XY]+);([0-9XY]+)\\)")[, 3],
      t_b1 = stringr::str_match(
        raw,
        "t\\([0-9XY]+;[0-9XY]+\\)\\(([^;)]+);([^)]+)\\)"
      )[, 2],
      t_b2 = stringr::str_match(
        raw,
        "t\\([0-9XY]+;[0-9XY]+\\)\\(([^;)]+);([^)]+)\\)"
      )[, 3]
    ) |>
    dplyr::filter(
      is_der,
      !is.na(t_c1),
      !is.na(der_inner),
      !stringr::str_detect(der_inner, ";")
    )

  if (nrow(der_t) == 0) {
    return(
      der_t |>
        dplyr::mutate(
          p1 = character(0),
          p2 = character(0),
          bb1 = character(0),
          bb2 = character(0),
          donor = character(0),
          balanced_pair = logical(0)
        )
    )
  }

  der_t |>
    dplyr::mutate(
      swap = .chrom_sort_key(t_c1) > .chrom_sort_key(t_c2),
      p1 = dplyr::if_else(swap, t_c2, t_c1),
      p2 = dplyr::if_else(swap, t_c1, t_c2),
      bb1 = dplyr::if_else(swap, t_b2, t_b1),
      bb2 = dplyr::if_else(swap, t_b1, t_b2),
      t_sig = paste0(
        "t(",
        p1,
        ";",
        p2,
        ")(",
        tidyr::replace_na(bb1, ""),
        ";",
        tidyr::replace_na(bb2, ""),
        ")"
      ),
      donor = der_inner
    ) |>
    dplyr::group_by(.pk_row_id, clone_id, t_sig) |>
    dplyr::mutate(
      balanced_pair = dplyr::if_else(
        p1 == p2,
        # Homologous t(a;a): both reciprocal products are der(a), so donor set
        # membership cannot distinguish one der from two. Require >= 2 ders for
        # the same signature before calling it a balanced reciprocal pair.
        dplyr::n() >= 2L,
        (p1 %in% donor) & (p2 %in% donor)
      )
    ) |>
    dplyr::ungroup()
}

# Classify each row as carrying balanced and/or unbalanced translocations.
#
# Balanced  := a bare t(...) token, OR a der()t() reciprocal pair where both
#              partner chromosomes appear as centromere donors within a clone.
# Unbalanced := a lone der(a)t(a;b) (the reciprocal der is absent), OR a
#              whole-arm der(a;b) (single derivative). A row may be both (mixed).
#
# These flags are independent of `derivative_chromosome`, which still fires.
classify_translocation_balance <- function(
  tokens_tbl,
  der_t = .der_translocations(tokens_tbl)
) {
  toks <- tokens_tbl |>
    dplyr::filter(!is.na(clone_id), aberr_raw != "") |>
    dplyr::mutate(
      raw = strip_bands(aberr_raw),
      is_der = stringr::str_detect(raw, "^\\+?i?der\\("),
      der_inner = stringr::str_match(raw, "^\\+?i?der\\(([^)]+)\\)")[, 2],
      is_bare_t = stringr::str_detect(raw, "^[?~]?t\\([0-9XY]")
    )

  empty <- tibble::tibble(
    .pk_row_id = integer(),
    balanced_translocation = integer(),
    unbalanced_translocation = integer()
  )

  bare_rows <- toks$.pk_row_id[toks$is_bare_t]

  # Whole-arm der(a;b): a single derivative, always unbalanced.
  whole_arm_rows <- toks$.pk_row_id[
    toks$is_der &
      !is.na(toks$der_inner) &
      stringr::str_detect(toks$der_inner, ";")
  ]

  der_balanced_rows <- der_t$.pk_row_id[der_t$balanced_pair]
  der_unbalanced_rows <- der_t$.pk_row_id[!der_t$balanced_pair]

  balanced_ids <- unique(c(bare_rows, der_balanced_rows))
  unbalanced_ids <- unique(c(der_unbalanced_rows, whole_arm_rows))
  all_ids <- sort(unique(c(balanced_ids, unbalanced_ids)))

  if (length(all_ids) == 0) {
    return(empty)
  }

  tibble::tibble(
    .pk_row_id = all_ids,
    balanced_translocation = as.integer(all_ids %in% balanced_ids),
    unbalanced_translocation = as.integer(all_ids %in% unbalanced_ids)
  )
}

# Derive implied partial losses from lone (unbalanced) der(a)t(a;b)(bp_a;bp_b):
#   - chromosome a loses material distal to bp_a -> arm(bp_a)
#   - partner b loses its centromere-side material  -> opposite arm of bp_b
# Emitted into dedicated `unbal_loss_<arm>` columns (one per chromosome arm,
# autosomes plus X/Y) plus a generic `unbal_partial_loss` flag. These are kept
# SEPARATE from del()/mono* signals -- an unbalanced-derived 5q loss is not
# conflated with a true del(5q).
#
# Scope: only SIMPLE single-junction ders are resolved -- exactly one two-partner
# t() whose centromere donor is one of the partners. Multi-junction chains
# (der with >1 t(), three-way t(a;b;c), or a der named after a non-partner
# chromosome) cannot be resolved by this two-partner model, so no loss is
# derived for them (they are still flagged unbalanced by the classifier).
# Copy-number-aware gains are out of scope (they need whole-karyotype reasoning).
derive_unbalanced_loss <- function(der_t) {
  cols <- paste0("unbal_loss_", .unbal_loss_arms)
  empty <- tibble::tibble(.pk_row_id = integer())
  for (nm in c(cols, "unbal_partial_loss")) {
    empty[[nm]] <- integer()
  }

  dt <- der_t |>
    dplyr::filter(
      !balanced_pair,
      stringr::str_count(raw, "t\\(") == 1L,
      donor == p1 | donor == p2
    )
  if (nrow(dt) == 0) {
    return(empty)
  }

  segs <- dt |>
    dplyr::mutate(
      donor_is_p1 = donor == p1,
      bp_a = dplyr::if_else(donor_is_p1, bb1, bb2),
      bp_b = dplyr::if_else(donor_is_p1, bb2, bb1),
      b_chrom = dplyr::if_else(donor_is_p1, p2, p1),
      a_seg = .arm_segment(donor, .arm_of(bp_a)),
      b_seg = .arm_segment(b_chrom, .opp_arm(.arm_of(bp_b)))
    )

  long <- dplyr::bind_rows(
    segs |> dplyr::transmute(.pk_row_id, seg = a_seg),
    segs |> dplyr::transmute(.pk_row_id, seg = b_seg)
  ) |>
    dplyr::filter(!is.na(seg))

  if (nrow(long) == 0) {
    return(empty)
  }

  any_loss <- long |>
    dplyr::distinct(.pk_row_id) |>
    dplyr::mutate(unbal_partial_loss = 1L)

  named <- long |>
    dplyr::filter(seg %in% .unbal_loss_arms) |>
    dplyr::distinct(.pk_row_id, seg) |>
    dplyr::mutate(col = paste0("unbal_loss_", seg), value = 1L) |>
    dplyr::select(.pk_row_id, col, value) |>
    tidyr::pivot_wider(
      names_from = col,
      values_from = value,
      values_fill = 0L
    )

  out <- any_loss |> dplyr::left_join(named, by = ".pk_row_id")
  for (nm in cols) {
    if (!nm %in% names(out)) out[[nm]] <- 0L
  }
  out |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(cols),
      ~ tidyr::replace_na(., 0L)
    ))
}

# Numeric sort key for a chromosome label (X -> 23, Y -> 24).
.chrom_sort_key <- function(chrom) {
  dplyr::case_when(
    chrom == "X" ~ 23L,
    chrom == "Y" ~ 24L,
    TRUE ~ suppressWarnings(as.integer(chrom))
  )
}

# First arm letter (p/q) of a band string, or NA.
.arm_of <- function(band) {
  stringr::str_extract(band, "^[pq]")
}

# Opposite chromosome arm.
.opp_arm <- function(arm) {
  dplyr::case_when(arm == "p" ~ "q", arm == "q" ~ "p", TRUE ~ NA_character_)
}

# Combine chromosome + arm into a segment key (e.g. "5q"), NA if arm unknown.
.arm_segment <- function(chrom, arm) {
  dplyr::if_else(is.na(arm) | is.na(chrom), NA_character_, paste0(chrom, arm))
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

# Provenance-tagged catalog of every parse_karyo() output column, in output
# order. This is the single source of truth from which parse_karyo() derives
# the full column list, the abnormality-flag subset, and character-column
# typing. Distinct from the *rules* schema (flag_name/regex/...) validated by
# validate_rules(): that describes the input rule table; this describes the
# output result table.
# Columns:
#   - column: output column name
#   - class:  provenance (meta, rule, aneuploidy, summary, balance, loss,
#             status). The `rule` rows depend on the active rule set.
#   - type:   storage type ("character" or "integer")
#   - description: one-line human description
.column_catalog <- function(rules) {
  chroms <- c(as.character(1:22), "X", "Y")
  row <- function(column, class, type, description) {
    tibble::tibble(
      column = column,
      class = class,
      type = type,
      description = description
    )
  }
  dplyr::bind_rows(
    row("original_karyotype", "meta", "character", "Raw input string"),
    row(
      "preprocessed_karyotype",
      "meta",
      "character",
      "Cleaned/normalized string that was parsed; NA for unfixable rows"
    ),
    row(
      "ploidy_category",
      "meta",
      "character",
      "Ploidy classification of the most abnormal clone"
    ),
    row(
      "chromosome_count",
      "meta",
      "integer",
      "Chromosome count from the most abnormal eligible clone"
    ),
    row(
      unique(rules$flag_name),
      "rule",
      "integer",
      "Aberration flag from the rule-matching pipeline"
    ),
    row(
      names(.general_flag_patterns),
      "general",
      "integer",
      "General structural-aberration flag; fires independently of specific rules"
    ),
    row(
      paste0("mono", chroms),
      "aneuploidy",
      "integer",
      "Monosomy flag derived from '-' tokens"
    ),
    row(
      paste0("tris", chroms),
      "aneuploidy",
      "integer",
      "Trisomy flag derived from '+' tokens"
    ),
    row(
      "normal_karyotype",
      "summary",
      "integer",
      "1 if 46,XX or 46,XY exactly"
    ),
    row(
      "comma_count_aberrations",
      "summary",
      "integer",
      "Aberration count (max across clones; idem-expanded)"
    ),
    row(
      "complex_karyotype",
      "summary",
      "integer",
      "1 if >= 3 distinct aberrations across clones"
    ),
    row(
      "monosomal_karyotype",
      "summary",
      "integer",
      "1 if monosomal-karyotype criteria are met"
    ),
    row(
      "mixed_ploidy",
      "summary",
      "integer",
      "1 if clones span different ploidy categories"
    ),
    row(
      "balanced_translocation",
      "balance",
      "integer",
      "1 if the row carries a balanced translocation"
    ),
    row(
      "unbalanced_translocation",
      "balance",
      "integer",
      "1 if the row carries an unbalanced translocation"
    ),
    row(
      paste0("unbal_loss_", .unbal_loss_arms),
      "loss",
      "integer",
      "Partial arm loss implied by an unbalanced der translocation"
    ),
    row(
      "unbal_partial_loss",
      "loss",
      "integer",
      "1 if any unbalanced-derived partial loss (OR of unbal_loss_*)"
    ),
    row(
      "total_metaphases",
      "meta",
      "integer",
      "Sum of bracket counts; NA if no brackets"
    ),
    row(
      "fixable_error",
      "status",
      "integer",
      "1 if the row had a fixable issue"
    ),
    row(
      "unfixable_error",
      "status",
      "integer",
      "1 if the row had an unfixable issue (row is all NA)"
    ),
    row(
      "chimeric_karyotype",
      "status",
      "integer",
      "1 if input contained a '//' chimeric separator"
    ),
    row(
      "chimeric_clone",
      "status",
      "character",
      "Which clone was parsed for a chimeric row (host/donor/NA)"
    )
  )
}

blank_rows <- function(original_karyotypes, all_output_cols, char_cols) {
  n <- length(original_karyotypes)
  out <- tibble::tibble(.rows = n)
  # Build in all_output_cols order so the column order matches a normal result.
  for (nm in all_output_cols) {
    out[[nm]] <- if (nm %in% char_cols) {
      rep(NA_character_, n)
    } else {
      rep(NA_integer_, n)
    }
  }
  out$original_karyotype <- original_karyotypes
  out
}
