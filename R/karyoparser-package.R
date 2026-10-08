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

.chrom_alt <- "X|Y|[1-9]|1[0-9]|2[0-2]"
# Any one- or two-digit number in a chromosome slot. Used by the dirty-pattern
# fixes and rule_match_text(), so a repair still applies to an impossible
# chromosome ('t(5:27)') and invalid_breakpoint then reports it.
.chrom_loose_alt <- "\\d{1,2}|X|Y"

# Whole-chromosome gain or loss: the entire token, optionally with a copy
# multiplier ('-Xx2'). Anchored at both ends so a leading copy count on an
# unidentified element ('+2mar', '+1~3r') is never read as a chromosome.
# Constitutional ('+21c') and uncertain ('+?8') forms deliberately don't
# match, so they set no mono/tris flag.
.aneuploidy_re <- paste0("^([+-])(", .chrom_alt, ")(?:x\\d+)?$")

# Optional sign, '?' and copy count (single or range) that can lead an
# unidentified ring, marker or double minute ('+1~3r', '+2mar', '-1mar').
# Prefix only: callers append the element and any copy multiplier ('x2').
.unidentified_count_re <- "^[+-]?\\??(?:\\d+(?:[~-]\\d+)?)?"
# The element itself, numbered ('mar1'), multiplied ('marx2') or uncertain.
.unidentified_element_re <- "(?:r|mar|dmin)\\d*(?:x\\d+)?\\??$"

# Every recognized shape of a '+'/'-' led token. Anything else is reported as
# unrecognized_gain_loss rather than silently ignored by the strict
# .aneuploidy_re.
.recognized_gain_loss_re <- c(
  aneuploidy = paste0("^\\??[+-]\\??(?:", .chrom_alt, ")(?:c|x\\d+)?\\??$"),
  unidentified = paste0(.unidentified_count_re, .unidentified_element_re),
  structural = "^\\??[+-]\\??[a-z][a-z ]*\\("
)
