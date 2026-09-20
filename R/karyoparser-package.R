#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs for dplyr/tidyr NSE column references
utils::globalVariables(c(
  ".pk_row_id",
  "aberr_norm",
  "aberr_raw",
  "aberrations",
  "a_seg",
  "autosomal_monosomies_sample",
  "b_chrom",
  "b_seg",
  "balanced_pair",
  "bb1",
  "bb2",
  "bp_a",
  "bp_b",
  "bracket",
  "chrom",
  "chromosome_count",
  "cleaned",
  "clone_id",
  "clone_str",
  "comma_count_aberrations",
  "content",
  "der_inner",
  "donor",
  "donor_is_p1",
  "first_is_sex",
  "first_token",
  "flag",
  "flag_name",
  "has_idem",
  "head",
  "issue_type",
  "is_composite",
  "is_cp",
  "is_der",
  "is_mono_tri",
  "is_range_karyotype",
  "fixable_error",
  "unfixable_error",
  "max_count",
  "monosomal_karyotype",
  "n",
  "n_partners",
  "n_unique_aberr",
  "normal_karyotype",
  "p1",
  "p2",
  "partners",
  "preprocessed_karyotype",
  "original_karyotype",
  "prefix",
  "raw",
  "raw_comma_count",
  "row_index",
  "seg",
  "stemline_aberr_count",
  "structural_aberrations_sample",
  "t_bands",
  "t_chroms",
  "t_sig",
  "tokens_nonempty",
  "total_metaphases",
  "value"
))


# Constants ----
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

# Structural aberration indicator prefixes, the single source of truth for
# these token strings. Used to ground auto-fixes that must tell a real
# aberration token from glued text (assess.R) and to build the `general_*`
# flags (flags.R, which loads after this file). All are paren-led except the
# bare ones below.
.aberr_indicators_paren <- c(
  psu_dic = "psu dic",
  idic = "idic",
  # Not a distinct general_* flag; flags.R folds it into general_derivative.
  ider = "ider",
  trp = "trp",
  dup = "dup",
  ins = "ins",
  inv = "inv",
  add = "add",
  del = "del",
  dic = "dic",
  der = "der",
  i = "i",
  r = "r",
  t = "t"
)
# Bare indicators: no breakpoint list, and both can carry a leading count
# ('2mar'), so their flag patterns must not require a left word boundary.
.aberr_indicators_bare <- c(mar = "mar", dmin = "dmin")
