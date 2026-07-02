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

# Extract karyotype and ID vectors from a data frame; shared by all three exported functions.
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
