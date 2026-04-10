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
  "karyo",
  "Karyo",
  "iscn",
  "ISCN"
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
  "id",
  "ID",
  "Id",
  "mrn",
  "MRN",
  "subject_id",
  "SubjectID",
  "subject"
)

# Extract karyotype and ID vectors from a plain data frame.
#
# Shared by check_karyo(), preprocess_karyo(), and parse_karyo() for the
# data-frame input path. Handles karyotype column detection (auto or explicit),
# type guard, and ID column detection (auto or explicit).
#
# @param x A data frame.
# @param karyotype_column Character, or NULL for auto-detection.
# @param id_column Character, or NULL for auto-detection.
# @param verbose Logical.
# @param caller Character string used in error messages.
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
    karyotype_col_name <- intersect(.karyotype_col_candidates, names(x))[1]
    if (is.na(karyotype_col_name)) {
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
    if (isTRUE(verbose)) {
      message("Auto-detected karyotype column: ", karyotype_col_name)
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
    id_col_name <- intersect(.id_col_candidates, names(x))[1]
    if (!is.na(id_col_name)) {
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
      id_col_name <- NULL
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
