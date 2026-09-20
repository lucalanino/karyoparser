# Columns are always declared, never inferred: guessing which column holds
# the karyotype risks silently parsing the wrong one when a data frame
# carries more than one plausible candidate (e.g. a raw `iscn` alongside a
# cleaned `karyotype`).

.quote_cols <- function(x) {
  if (length(x) == 0L) {
    return("none")
  }
  paste(sprintf("\"%s\"", x), collapse = ", ")
}

# Duplicate ids are a downstream problem, not a karyoparser one: rows are
# matched on the karyotype string, never on the id, so parsing is unaffected.
# Warn rather than stop -- serial karyotypes for one patient (diagnosis,
# post-induction, relapse) are a legitimate shape.
.warn_duplicate_ids <- function(id_values, id_col_name) {
  dup_mask <- duplicated(id_values) | duplicated(id_values, fromLast = TRUE)
  if (!any(dup_mask)) {
    return(invisible(NULL))
  }
  dup_values <- unique(id_values[dup_mask])
  n_na <- sum(is.na(dup_values))
  shown <- utils::head(dup_values[!is.na(dup_values)], 5L)
  detail <- if (length(shown) > 0L) {
    sprintf(
      " (%s%s)",
      .quote_cols(shown),
      if (length(dup_values) - n_na > length(shown)) ", ..." else ""
    )
  } else {
    ""
  }
  warning(
    sprintf(
      paste0(
        "ID column '%s' has %d duplicated value(s) across %d row(s)%s.%s",
        "\nkaryoparser is unaffected -- rows are matched on the karyotype",
        " string, not the id -- but joins on this column downstream may fan",
        " out or fail."
      ),
      id_col_name,
      length(dup_values),
      sum(dup_mask),
      detail,
      if (n_na > 0L) " Duplicates include NA." else ""
    ),
    call. = FALSE
  )
  invisible(NULL)
}

.extract_df_input <- function(
  x,
  karyotype_column,
  id_column,
  verbose,
  caller
) {
  if (is.null(karyotype_column)) {
    stop(
      sprintf(
        paste(
          "`karyotype_column` must be given when %s() is passed a data frame.",
          "\nAvailable columns: %s.",
          "\ne.g. %s(data, karyotype_column = \"<col>\")"
        ),
        caller,
        .quote_cols(names(x)),
        caller
      ),
      call. = FALSE
    )
  }
  if (!karyotype_column %in% names(x)) {
    stop(
      sprintf(
        "Karyotype column '%s' not found in input. Available columns: %s.",
        karyotype_column,
        .quote_cols(names(x))
      ),
      call. = FALSE
    )
  }
  karyotype_col_name <- karyotype_column
  if (isTRUE(verbose)) {
    message("Using karyotype column: ", karyotype_col_name)
  }

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

  id_col_name <- NULL
  if (!is.null(id_column)) {
    if (!id_column %in% names(x)) {
      stop(
        sprintf(
          "ID column '%s' not found in input. Available columns: %s.",
          id_column,
          .quote_cols(names(x))
        ),
        call. = FALSE
      )
    }
    id_col_name <- id_column
    if (isTRUE(verbose)) {
      message("Using ID column: ", id_col_name)
    }
  } else if (isTRUE(verbose)) {
    message("No ID column given. Use id_column= to carry one through.")
  }

  id_values <- if (!is.null(id_col_name)) {
    as.character(x[[id_col_name]])
  } else {
    NULL
  }
  if (!is.null(id_values)) {
    .warn_duplicate_ids(id_values, id_col_name)
  }

  list(
    raw_vec = raw_vec,
    karyotype_col_name = karyotype_col_name,
    id_values = id_values,
    id_col_name = id_col_name
  )
}
