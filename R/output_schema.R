# Single source of truth for parse_karyo()'s output column list, order, and typing.
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
      "chromosome_count",
      "meta",
      "integer",
      "Chromosome count from the most abnormal eligible clone"
    ),
    row(
      .sanitize_flag_name(unique(rules$flag_name)),
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
      paste0("unbal_partial_loss_", .unbal_partial_loss_arms),
      "loss",
      "integer",
      "Partial arm loss implied by an unbalanced der translocation"
    ),
    row(
      "unbal_partial_loss",
      "loss",
      "integer",
      "1 if any unbalanced-derived partial loss (OR of unbal_partial_loss_*)"
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
