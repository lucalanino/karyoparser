#' Parse ISCN Karyotype Strings
#'
#' Parses ISCN karyotype notation into structured binary features for analysis.
#' Extracts specific translocations, deletions, monosomies, trisomies, and other
#' chromosomal aberrations according to configurable rules. Performs no
#' cleaning or normalization itself -- that is `preprocess_karyo()`'s job --
#' but `on_issues = "preprocess"` invokes the same checking/fixing logic as
#' `check_karyo()`/`preprocess_karyo()` internally, so a single call is
#' usually enough without chaining all three explicitly.
#'
#' @param karyotypes A character vector of karyotype strings, a data frame
#'   containing a karyotype column, or a `karyo_preprocessed` tibble returned
#'   by `preprocess_karyo()`. When a `karyo_preprocessed` object is passed, the
#'   `preprocessed` column is used automatically and cached issue indices and id
#'   column are reused -- no additional arguments required.
#' @param rules A `karyo_rules` object (default: [myeloid_rules]), or a list
#'   of `karyo_rules` objects to combine (e.g. `list(myeloid_rules,
#'   lymphoid_rules)`). Rules already fire independently within a single
#'   table -- multiple rules can match the same token -- so combining tables
#'   just widens the pool. `flag_name` must be unique across the combined
#'   tables: [validate_rules()] errors if the same `flag_name` appears more
#'   than once, since there is no cross-table priority to resolve the
#'   conflict. Use [validate_rules()] to validate and convert a custom data
#'   frame into an accepted rules object. Each rule's output column is a
#'   sanitized version of its `flag_name` (non-alphanumeric characters
#'   replaced by `_`), e.g. `"t(9;22)(q34;q11)"` becomes `t_9_22_q34_q11`.
#' @param karyotype_column Character. Name of the karyotype column when input is
#'   a plain data frame. If `NULL` (default), auto-detected from common names
#'   (`karyotype`, `iscn`, etc.). Ignored for character vector or
#'   `karyo_preprocessed` input.
#' @param id_column Character. Name of the id column when input is a data frame.
#'   If `NULL` (default), auto-detected from common names (`sample_id`,
#'   `patient_id`, `id`, `mrn`, etc.) and used silently. Detection follows the
#'   candidate list's own order rather than the column order, so `sample_id`
#'   wins over `mrn` whatever their positions. The id column is placed first
#'   in the output.
#'
#'   For `karyo_preprocessed` input the id is inherited automatically and
#'   `id_column` is rarely needed. `preprocess_karyo()` returns a narrow
#'   tibble -- its id column plus `original`, `preprocessed`, and `status` --
#'   so every other column of the source data frame is already gone, and
#'   naming one here is an error rather than an override. To key the output
#'   on a different column, name it upstream (`check_karyo()` or
#'   `preprocess_karyo(..., id_column = )`), or attach it to the preprocessed
#'   tibble yourself before parsing. The three structural column names are
#'   reserved and rejected as `id_column`.
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
#' @param columns Character vector selecting which output-column classes to
#'   include, or `NULL` (default) for the full output. `original_karyotype`,
#'   `preprocessed_karyotype`, `fixable_error`, `unfixable_error`,
#'   `chimeric_karyotype`, and `chimeric_clone` are always included regardless
#'   of `columns`. Valid values (any combination):
#'   - `"classification"`: `normal_karyotype`, `complex_karyotype`,
#'     `monosomal_karyotype`
#'   - `"rule"`: rule-matching aberration flags (from `rules`)
#'   - `"general"`: disease-agnostic structural flags (`general_*`,
#'     `balanced_translocation`, `unbalanced_translocation`)
#'   - `"aneuploidy"`: `mono*`/`tris*` monosomy/trisomy flags
#'   - `"summary"`: `comma_count_aberrations`, `chromosome_count`,
#'     `total_metaphases`
#'   - `"der_loss"`: `unbal_partial_loss_<arm>` and `unbal_partial_loss`
#' @param min_metaphases Numeric. Minimum metaphase count for a clone to be
#'   eligible to supply `chromosome_count`. Default `2`. Clones whose bracket
#'   reports fewer metaphases are skipped, so a small side clone does not
#'   outvote the main one; if no clone in a row clears the threshold,
#'   `chromosome_count` is `NA` for that row. The threshold applies only to
#'   rows where at least one clone carries a metaphase bracket -- a row with
#'   no brackets anywhere uses all of its clones regardless. Raise it to
#'   ignore more small clones (`5` is a common cytogenetics convention), or
#'   set `0` to consider every clone that has a count.
#'
#' @return A tibble with columns:
#'   - original_karyotype: Input karyotype string (always the raw input)
#'   - preprocessed_karyotype: `normalize_iscn()` output of the
#'     karyotype string, consistent across all `on_issues` modes. In
#'     `"preprocess"` mode this is also dirty-fixed; in `"warn"`/`"stop"` modes
#'     only `normalize_iscn()` is applied (chimeric rows are clone-selected in
#'     every mode). `NA_character_` for issue rows.
#'   - chromosome_count: Integer chromosome count, taken from the most
#'     abnormal eligible clone in the row (the one whose count is furthest
#'     from 46). Eligibility is governed by `min_metaphases`; `NA` when no
#'     clone in the row qualifies, which can happen even though
#'     `total_metaphases` is non-`NA`.
#'   - One column per aberration flag (0/1 binary); rule-based columns are
#'     named after a sanitized `flag_name` (see `rules` above)
#'   - monoX, monoY, mono1-22: Monosomy flags for each chromosome
#'   - trisX, trisY, tris1-22: Trisomy flags for each chromosome
#'   - normal_karyotype: 1 if 46,XX or 46,XY, else 0
#'   - total_metaphases: Count from bracket notation
#'   - comma_count_aberrations: Number of comma-separated aberrations
#'   - distinct_aberrations: Number of distinct normalized aberration
#'     tokens pooled across all clones, excluding sex complements,
#'     `idem`/`sl`, and range markers. This is the count
#'     `complex_karyotype` thresholds at 3.
#'   - complex_karyotype: 1 if >=3 unique aberrations, else 0
#'   - monosomal_karyotype: 1 if the row has two or more autosomal
#'     monosomies, or one autosomal monosomy plus at least one structural
#'     aberration (any `general_*` flag except `general_marker` -- a lone
#'     marker does not count). Two carve-outs: sex-chromosome monosomies
#'     (`monoX`/`monoY`) never count toward either criterion, and a
#'     CBF-AML row (`t(8;21)(q22;q22)`, `inv(16)(p13q22)`, or
#'     `t(16;16)(p13;q22)`) is forced to 0 per clinical guidelines even
#'     when the criteria are met.
#'   - balanced_translocation: 1 if the row carries a balanced translocation --
#'     a bare `t(...)` token, or a reciprocal der pair where both partner
#'     chromosomes appear as centromere donors (e.g.
#'     `der(5)t(5;17)...,der(17)t(5;17)...`). Independent of
#'     `general_derivative`. A row may be both balanced and unbalanced.
#'   - unbalanced_translocation: 1 if the row carries an unbalanced
#'     translocation -- a lone der whose reciprocal set is incomplete (e.g.
#'     `der(a)t(a;b)` or a lone `der(a)t(a;b;c)` three-way der) or a whole-arm
#'     `der(a;b)`. The specific fusion flag still fires (e.g. `der(9)t(9;22)`
#'     keeps `t_9_22_q34_q11 = 1`). A der is balanced instead only when the
#'     full reciprocal set is present (every partner appears as a der).
#'   - unbal_partial_loss_<arm>: one column per chromosome arm
#'     (`unbal_partial_loss_1p`, `unbal_partial_loss_1q`, ...,
#'     `unbal_partial_loss_22q`, `unbal_partial_loss_Xp`, `unbal_partial_loss_Xq`,
#'     `unbal_partial_loss_Yp`, `unbal_partial_loss_Yq`), set to 1 when an
#'     unbalanced der translocation implies partial loss of that arm. Kept
#'     SEPARATE from `del(...)`/`mono*` -- an unbalanced-derived 5q loss does not
#'     set `del_5q`. Derivation needs explicit breakpoints and is limited to
#'     simple single-junction two-partner `der(a)t(a;b)`; multi-junction chains
#'     and three-way `t(a;b;c)` derivatives are flagged unbalanced but derive no
#'     loss. Copy-number-aware gains are out of scope.
#'   - unbal_partial_loss: 1 if an unbalanced der translocation implies any
#'     partial loss (OR across all `unbal_partial_loss_*` columns).
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
#' @seealso [check_karyo()], [preprocess_karyo()], `vignette("data-cleaning")`
#' @export
parse_karyo <- function(
  karyotypes,
  rules = myeloid_rules,
  karyotype_column = NULL,
  id_column = NULL,
  verbose = FALSE,
  on_issues = c("stop", "preprocess", "warn"),
  on_chimeric = c("default", "host", "donor"),
  columns = NULL,
  min_metaphases = 2
) {
  on_issues <- match.arg(on_issues)
  on_chimeric_provided <- !missing(on_chimeric)
  on_chimeric <- match.arg(on_chimeric)

  if (
    !is.numeric(min_metaphases) ||
      length(min_metaphases) != 1L ||
      is.na(min_metaphases) ||
      min_metaphases < 0
  ) {
    stop(
      "`min_metaphases` must be a single non-negative number.",
      call. = FALSE
    )
  }

  # Validate rules early (needed for empty result structure) ------------------
  # `rules` may be a single karyo_rules object, or a list of karyo_rules
  # objects to combine (e.g. list(myeloid_rules, lymphoid_rules)). Excluding
  # data frames here matters because a data.frame is itself a list -- without
  # it, a plain (unvalidated) rules data frame would be mistaken for a list of
  # rules tables instead of hitting the karyo_rules class check below.
  if (!is.data.frame(rules) && is.list(rules)) {
    if (length(rules) == 0L) {
      stop(
        "`rules` must contain at least one `karyo_rules` object.",
        call. = FALSE
      )
    }
    is_karyo_rules <- vapply(rules, inherits, logical(1), what = "karyo_rules")
    if (!all(is_karyo_rules)) {
      stop(
        "Every element of `rules` must be a `karyo_rules` object. ",
        "Use `validate_rules()` to validate a custom rules table.",
        call. = FALSE
      )
    }
    rules <- validate_rules(dplyr::bind_rows(lapply(rules, tibble::as_tibble)))
  }
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

  # ---- Output-column filtering (width knob) ----------------------------------
  # `columns` selects which output-column classes to keep, using the `class`
  # field from .column_catalog(). "meta" and "status" (row identity plus
  # error/chimeric bookkeeping) are always kept regardless of `columns` --
  # every caller needs them to identify rows and interpret the rest, mirroring
  # how `id_column` is always placed first regardless of other arguments.
  filterable_classes <- setdiff(unique(catalog$class), c("meta", "status"))
  if (!is.null(columns)) {
    unknown_classes <- setdiff(columns, filterable_classes)
    if (length(unknown_classes) > 0) {
      stop(
        sprintf(
          "Unknown `columns` value(s): %s. Valid classes: %s.",
          paste(unknown_classes, collapse = ", "),
          paste(filterable_classes, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    catalog <- catalog[
      catalog$class %in% c("meta", "status", columns),
      ,
      drop = FALSE
    ]
  }

  all_output_cols <- catalog$column
  # chromosome_count/total_metaphases excluded: NA is meaningful for them, not backfillable to 0.
  abnormality_names <- catalog$column[
    catalog$class %in%
      c("rule", "general", "aneuploidy", "classification", "der_loss") |
      catalog$column %in% c("comma_count_aberrations", "distinct_aberrations")
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
      # preprocess_karyo() returns a narrow tibble -- its id column plus
      # original/preprocessed/status -- so any other column of the source
      # data frame is already gone by the time we get here, and the only
      # real fix is to name it upstream. The three structural columns are
      # reserved: accepting one would silently fill the output's id column
      # with karyotype strings or cleaning statuses.
      reserved <- c("original", "preprocessed", "status")
      if (id_column %in% reserved) {
        stop(
          sprintf(
            paste(
              "`id_column` cannot be \"%s\": that is a structural column of",
              "karyo_preprocessed, not an identifier.",
              "\nUsing it would fill the id column with karyotype strings or",
              "cleaning statuses."
            ),
            id_column
          ),
          call. = FALSE
        )
      }
      if (!id_column %in% names(karyotypes)) {
        available <- setdiff(names(karyotypes), reserved)
        stop(
          sprintf(
            paste(
              "ID column \"%s\" not found in karyo_preprocessed input.",
              "Available: %s.",
              "\npreprocess_karyo() keeps only its id column plus",
              "`original`/`preprocessed`/`status`, so a column dropped",
              "upstream cannot be named here.",
              "\nRe-run check_karyo() or preprocess_karyo() with",
              "id_column = \"%s\" to carry it through, or attach the column",
              "to the preprocessed tibble before parsing."
            ),
            id_column,
            if (length(available) > 0) {
              paste(sprintf("\"%s\"", available), collapse = ", ")
            } else {
              "none"
            },
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
        # Blank unfixable rows plus any chimeric row whose selected clone is NA
        # (e.g. a zero-host chimera under on_chimeric = "host"): output is NA,
        # but such policy-NA rows stay classified as fixable, not unfixable.
        chimeric_na_rows <- intersect(which(is.na(proc)), chimeric_all_indices)
        issue_row_indices <- unique(c(unfixable_row_indices, chimeric_na_rows))

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
  sample_meta <- build_sample_meta(dedup_df, min_metaphases = min_metaphases)
  tokens <- build_clone_tokens(sample_meta)
  flags <- match_rules(tokens, rules, rule_flag_names)
  aneuploidy <- compute_aneuploidy(tokens, chroms)
  general_flags <- compute_general_flags(tokens)
  comma_counts <- compute_comma_counts(tokens)
  unique_aberr <- compute_unique_counts(tokens)
  der_translocations <- .der_translocations(tokens)
  balance_flags <- classify_translocation_balance(tokens, der_translocations)
  unbal_partial_loss <- derive_unbalanced_loss(der_translocations)

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
  cbf_flags <- .sanitize_flag_name(
    c("t(8;21)(q22;q22)", "inv(16)(p13q22)", "t(16;16)(p13;q22)")
  )
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
          distinct_aberrations = as.integer(n_unique_aberr),
          complex_karyotype = as.integer(n_unique_aberr >= 3)
        ),
      by = ".pk_row_id"
    ) |>
    dplyr::left_join(comma_counts, by = ".pk_row_id") |>
    dplyr::left_join(balance_flags, by = ".pk_row_id") |>
    dplyr::left_join(unbal_partial_loss, by = ".pk_row_id") |>
    dplyr::left_join(
      sample_meta |>
        dplyr::select(
          .pk_row_id,
          normal_karyotype,
          chromosome_count,
          total_metaphases
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
      )
    ) |>
    dplyr::left_join(dedup_df, by = ".pk_row_id") |>
    dplyr::rename(preprocessed_karyotype = original_karyotype) |>
    dplyr::relocate(preprocessed_karyotype, .before = 1)

  keep_cols <- intersect(all_output_cols, names(result))

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
