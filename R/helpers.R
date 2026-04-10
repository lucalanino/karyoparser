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

strip_bands <- function(band_str) {
  stringr::str_replace_all(band_str, "([pq]\\d+)\\.\\d+", "\\1")
}

normalize_token <- function(x) {
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "(^\\?|\\?$|~)", "")
  x <- strip_bands(x)
  out <- x
  # vectorized classification
  is_marker <- tidyr::replace_na(stringr::str_detect(x, "^[+-]?\\d*mar"), FALSE)
  out[is_marker] <- "marker_chromosomes"

  is_der <- tidyr::replace_na(stringr::str_detect(x, "^\\+?i?der\\("), FALSE)
  if (any(is_der)) {
    out[is_der] <- paste0(
      "der_",
      stringr::str_extract(x[is_der], "i?der\\([^)]+\\)")
    )
  }

  is_dic <- tidyr::replace_na(
    stringr::str_detect(x, "^(dic|idic|psu dic)\\("),
    FALSE
  )
  if (any(is_dic)) {
    head <- stringr::str_extract(x[is_dic], "^(dic|idic|psu dic)")
    paren <- stringr::str_extract(x[is_dic], "\\([^)]+\\)")
    out[is_dic] <- paste0(head, "_", paren)
  }
  out
}

# Candidate column names for auto-detection of the karyotype column in data
# frames passed to check_karyo(), preprocess_karyo(), or parse_karyo().
.karyotype_col_candidates <- c(
  "karyotype",
  "Karyotype",
  "KARYOTYPE",
  "karyotype_string",
  "karyo",
  "Karyo",
  "iscn",
  "ISCN",
  "iscn_string"
)

# Candidate column names for auto-detection of the id column in data frames.
.id_col_candidates <- c(
  "sample_id",
  "Sample_ID",
  "SampleID",
  "sample",
  "patient_id",
  "Patient_ID",
  "PatientID",
  "patient",
  "case_id",
  "CaseID",
  "case",
  "id",
  "ID",
  "Id",
  "mrn",
  "MRN",
  "accession_id",
  "accession",
  "specimen_id",
  "subject_id",
  "SubjectID",
  "subject"
)

# Abort helper: errors with a copy-pasteable manual-call hint.
# @keywords internal
.abort_column_selection <- function(caller) {
  stop(
    sprintf(
      paste0(
        "Column selection aborted.\n",
        "Specify both columns explicitly when calling `%s()`:\n",
        "  %s(df, karyotype_column = \"<col>\", id_column = \"<col>\")\n",
        "  # use id_column = NULL to omit the ID column"
      ),
      caller,
      caller
    ),
    call. = FALSE
  )
}

# Interactive column confirmation helper.
#
# Shows the auto-detected column (with a short preview) and prompts the user to
# accept, pick a different column, or abort. Only called when interactive() is
# TRUE.
#
# @param x        Data frame.
# @param detected Auto-detected column name (character), or NA_character_ when
#                 none was found.
# @param role     "karyotype" or "ID" — used in prompt text.
# @param caller   Calling function name — used in the abort hint.
# @param allow_none If TRUE, offer a "none" option to skip the column (for the
#                 ID column).
# @return Chosen column name (character), or NULL when user picks "none".
# @keywords internal
.confirm_column <- function(x, detected, role, caller, allow_none = FALSE) {
  cols <- names(x)

  if (!is.na(detected)) {
    # Auto-detected — ask user to confirm or redirect
    preview <- paste(
      sprintf('"%s"', utils::head(as.character(x[[detected]]), 3)),
      collapse = ", "
    )
    opts <- if (allow_none) {
      "Enter/y = use,  n = pick another,  none = skip,  stop = abort"
    } else {
      "Enter/y = use,  n = pick another,  stop = abort"
    }
    message(sprintf(
      '\nAuto-detected %s column: "%s"\n  Preview: %s\n  %s',
      role,
      detected,
      preview,
      opts
    ))
    answer <- trimws(readline("> "))

    if (answer == "" || tolower(answer) %in% c("y", "yes")) {
      return(detected)
    }
    if (allow_none && tolower(answer) == "none") {
      return(NULL)
    }
    if (tolower(answer) == "stop") {
      .abort_column_selection(caller)
    }
    # If user typed a column name directly (not "n"/"no"), try to use it
    if (!tolower(answer) %in% c("n", "no")) {
      if (answer %in% cols) {
        return(answer)
      }
      stop(sprintf("Column '%s' not found in input.", answer), call. = FALSE)
    }
    # User said "n" — fall through to pick-from-list
    candidate_cols <- setdiff(cols, detected)
  } else if (allow_none) {
    # ID column not found — offer to pick one or skip
    message(sprintf("\nNo %s column auto-detected.", role))
    candidate_cols <- cols
  } else {
    # Karyotype column not found — must pick one
    message(sprintf(
      '\nCould not auto-detect %s column in `%s()`.',
      role,
      caller
    ))
    candidate_cols <- cols
  }

  # Pick-from-list
  message("  Available columns: ", paste(candidate_cols, collapse = ", "))
  opts2 <- if (allow_none) {
    "Enter column name,  none = skip,  stop = abort"
  } else {
    "Enter column name,  stop = abort"
  }
  message("  ", opts2)
  answer2 <- trimws(readline("> "))

  if (tolower(answer2) == "stop") {
    .abort_column_selection(caller)
  }
  if (allow_none && (tolower(answer2) == "none" || answer2 == "")) {
    return(NULL)
  }
  if (!allow_none && answer2 == "") {
    .abort_column_selection(caller)
  }
  if (answer2 %in% cols) {
    return(answer2)
  }
  stop(sprintf("Column '%s' not found in input.", answer2), call. = FALSE)
}

# Extract karyotype and ID vectors from a plain data frame.
#
# Shared by check_karyo(), preprocess_karyo(), and parse_karyo() for the
# data-frame input path. Handles karyotype column detection (auto or explicit),
# type guard, and ID column detection (auto or explicit).
#
# When neither column argument is supplied and the session is interactive,
# auto-detected columns are confirmed with the user via readline() before use.
# In non-interactive sessions (scripts, CI), auto-detected columns are used
# silently (verbose=TRUE prints a notice).
#
# @param x A data frame.
# @param karyotype_column Character, or NULL for auto-detection.
# @param id_column Character, or NULL for auto-detection.
# @param verbose Logical.
# @param caller Character string used in error/prompt messages.
# @return Named list: raw_vec, karyotype_col_name, id_values, id_col_name.
# @keywords internal
.extract_df_input <- function(
  x,
  karyotype_column,
  id_column,
  verbose,
  caller
) {
  # --- Karyotype column -------------------------------------------------------
  if (is.null(karyotype_column)) {
    detected_karyo <- intersect(.karyotype_col_candidates, names(x))[1]
    if (interactive()) {
      karyotype_col_name <- .confirm_column(
        x,
        detected_karyo,
        "karyotype",
        caller,
        allow_none = FALSE
      )
    } else {
      if (is.na(detected_karyo)) {
        stop(
          sprintf(
            "Could not auto-detect karyotype column in `%s()`. ",
            caller
          ),
          "Tried: ",
          paste(.karyotype_col_candidates, collapse = ", "),
          ". Pass `karyotype_column = \"<col>\"` explicitly.",
          call. = FALSE
        )
      }
      karyotype_col_name <- detected_karyo
      if (isTRUE(verbose)) {
        message("Auto-detected karyotype column: ", karyotype_col_name)
      }
    }
  } else {
    if (!karyotype_column %in% names(x)) {
      stop(
        sprintf("Column '%s' not found in input.", karyotype_column),
        call. = FALSE
      )
    }
    karyotype_col_name <- karyotype_column
    if (isTRUE(verbose)) {
      message("Using karyotype column: ", karyotype_col_name)
    }
  }

  # Type guard: karyotype column must be character (or factor, which coerces)
  col_val <- x[[karyotype_col_name]]
  if (!is.character(col_val)) {
    warning(
      sprintf(
        "Column '%s' is %s, not character \u2014 coercing with as.character(). ",
        karyotype_col_name,
        class(col_val)[1]
      ),
      "Verify this column contains ISCN karyotype strings.",
      call. = FALSE
    )
  }
  raw_vec <- as.character(col_val)

  # --- ID column --------------------------------------------------------------
  if (is.null(id_column)) {
    detected_id <- intersect(.id_col_candidates, names(x))[1]
    if (interactive()) {
      id_col_name <- .confirm_column(
        x,
        detected_id,
        "ID",
        caller,
        allow_none = TRUE
      )
    } else {
      id_col_name <- if (!is.na(detected_id)) detected_id else NULL
      if (!is.null(id_col_name)) {
        if (isTRUE(verbose)) {
          message(
            "Auto-detected ID column: ",
            id_col_name,
            ". Use id_column= to specify explicitly."
          )
        }
      } else {
        if (isTRUE(verbose)) {
          message("No ID column detected. Use id_column= to specify one.")
        }
      }
    }
  } else {
    if (!id_column %in% names(x)) {
      stop(
        sprintf("ID column '%s' not found in input.", id_column),
        call. = FALSE
      )
    }
    id_col_name <- id_column
    if (isTRUE(verbose)) {
      message("Using ID column: ", id_col_name)
    }
  }

  id_values <- if (!is.null(id_col_name)) {
    as.character(x[[id_col_name]])
  } else {
    NULL
  }

  list(
    raw_vec = raw_vec,
    karyotype_col_name = karyotype_col_name,
    id_values = id_values,
    id_col_name = id_col_name
  )
}
