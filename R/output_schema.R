# Single source of truth for parse_karyo()'s output column list, order, and typing.
.column_catalog <- function(rules) {
  chroms <- c(as.character(1:22), "X", "Y")
  row <- function(column, class, type) {
    tibble::tibble(column = column, class = class, type = type)
  }
  dplyr::bind_rows(
    row("original_karyotype", "meta", "character"),
    row("preprocessed_karyotype", "meta", "character"),
    row("normal_karyotype", "classification", "integer"),
    row("complex_karyotype", "classification", "integer"),
    row("monosomal_karyotype", "classification", "integer"),
    row(.sanitize_flag_name(unique(rules$flag_name)), "rule", "integer"),
    row(names(.general_flag_patterns), "general", "integer"),
    row("balanced_translocation", "general", "integer"),
    row("unbalanced_translocation", "general", "integer"),
    row(paste0("mono", chroms), "aneuploidy", "integer"),
    row(paste0("tris", chroms), "aneuploidy", "integer"),
    row("comma_count_aberrations", "summary", "integer"),
    row("chromosome_count", "summary", "integer"),
    row("total_metaphases", "summary", "integer"),
    row(
      paste0("unbal_partial_loss_", .unbal_partial_loss_arms),
      "der_loss",
      "integer"
    ),
    row("unbal_partial_loss", "der_loss", "integer"),
    row("fixable_error", "status", "integer"),
    row("unfixable_error", "status", "integer"),
    row("chimeric_karyotype", "status", "integer"),
    row("chimeric_clone", "status", "character")
  )
}

blank_rows <- function(original_karyotypes, all_output_cols, char_cols) {
  n <- length(original_karyotypes)
  out <- tibble::tibble(.rows = n)
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
