#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs for dplyr/tidyr NSE column references
utils::globalVariables(c(
  ".pk_row_id",
  "aberr_norm",
  "aberr_raw",
  "aberrations",
  "autosomal_monosomies_sample",
  "bracket",
  "chrom",
  "chromosome_count",
  "cleaned",
  "clone_id",
  "clone_str",
  "comma_count_aberrations",
  "content",
  "counts_for_monosomal",
  "first_is_sex",
  "first_token",
  "flag",
  "flag_name",
  "has_idem",
  "head",
  "idem_invalid",
  "issue_type",
  "is_composite",
  "is_cp",
  "is_mono_tri",
  "is_range_karyotype",
  "fixable_error",
  "unfixable_error",
  "max_count",
  "mixed_ploidy",
  "monosomal_karyotype",
  "n",
  "n_unique_aberr",
  "normal_karyotype",
  "preprocessed_karyotype",
  "original_karyotype",
  "ploidy_result",
  "prefix",
  "priority",
  "raw_comma_count",
  "row_index",
  "stemline_aberr_count",
  "structural_aberrations_sample",
  "tokens_nonempty",
  "total_metaphases",
  "value"
))


# ---- Constants ---------------------------------------------------------------
.karyoparser_version <- as.character(utils::packageVersion("karyoparser"))

.sex_complements <- c(
  "X",
  "Y",
  "XX",
  "XY",
  "XXX",
  "XXY",
  "XYY",
  "XXYY",
  "XXXY",
  "XXXX",
  "XXXXY"
)
