# Parse ISCN Karyotype Strings

Parses ISCN karyotype notation into structured binary features for
analysis. Extracts specific translocations, deletions, monosomies,
trisomies, and other chromosomal aberrations according to configurable
rules. Performs no cleaning or normalization itself – that is
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)'s
job – but `on_issues = "preprocess"` invokes the same checking/fixing
logic as
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md)/[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
internally, so a single call is usually enough without chaining all
three explicitly.

## Usage

``` r
parse_karyo(
  karyotypes,
  rules = myeloid_rules,
  karyotype_column = NULL,
  id_column = NULL,
  verbose = FALSE,
  on_issues = c("stop", "preprocess", "warn"),
  on_chimeric = c("default", "host", "donor"),
  columns = NULL,
  min_metaphases = 2
)
```

## Arguments

- karyotypes:

  A character vector of karyotype strings, a data frame containing a
  karyotype column, or a `karyo_preprocessed` tibble returned by
  [`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md).
  When a `karyo_preprocessed` object is passed, the `preprocessed`
  column is used automatically and cached issue indices and id column
  are reused – no additional arguments required.

- rules:

  A `karyo_rules` object (default:
  [myeloid_rules](https://lucalanino.github.io/karyoparser/reference/myeloid_rules.md)),
  or a list of `karyo_rules` objects to combine (e.g.
  `list(myeloid_rules, lymphoid_rules)`). Rules already fire
  independently within a single table – multiple rules can match the
  same token – so combining tables just widens the pool. `flag_name`
  must be unique across the combined tables:
  [`validate_rules()`](https://lucalanino.github.io/karyoparser/reference/validate_rules.md)
  errors if the same `flag_name` appears more than once, since there is
  no cross-table priority to resolve the conflict. Use
  [`validate_rules()`](https://lucalanino.github.io/karyoparser/reference/validate_rules.md)
  to validate and convert a custom data frame into an accepted rules
  object. Each rule's output column is a sanitized version of its
  `flag_name` (non-alphanumeric characters replaced by `_`), e.g.
  `"t(9;22)(q34;q11)"` becomes `t_9_22_q34_q11`.

- karyotype_column:

  Character. Name of the karyotype column. **Required** when input is a
  plain data frame – columns are never guessed, since a frame carrying
  more than one plausible candidate (a raw `iscn` beside a cleaned
  `karyotype`) would otherwise be parsed from the wrong one silently.
  Ignored for character vector or `karyo_preprocessed` input.

- id_column:

  Character. Name of the id column. Optional: `NULL` (default) means the
  data has no identifier, and one is never inferred. When given, the
  column is placed first in the output and carried through the pipeline.
  Duplicate ids are warned about, not rejected: rows are matched on the
  karyotype string and never on the id, so parsing is unaffected – but
  downstream joins on the column may fan out.

  For `karyo_preprocessed` input the id is inherited automatically and
  `id_column` is rarely needed.
  [`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
  returns a narrow tibble – its id column plus `original`,
  `preprocessed`, and `status` – so every other column of the source
  data frame is already gone, and naming one here is an error rather
  than an override. To key the output on a different column, name it
  upstream
  ([`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md)
  or `preprocess_karyo(..., id_column = )`), or attach it to the
  preprocessed tibble yourself before parsing. The three structural
  column names are reserved and rejected as `id_column`.

- verbose:

  Logical. If `TRUE`, prints column detection and parsing summary
  messages. Default `FALSE`.

- on_issues:

  Character string specifying how to handle **dirty markers**
  (formatting artifacts such as leading dots, HTML entities, missing sex
  comma) and **structural errors** (e.g. no chromosome count,
  `Updated ISCN` marker). Chimeric handling is independent of this
  argument – see `on_chimeric`. Structural errors are always unfixable
  and return NA regardless of `on_issues`.

  Values:

  - `"stop"` (default): Raise an error immediately if any non-chimeric
    issue is found, listing fixable and unfixable issue types and
    counts, and suggesting next steps. No fixing is attempted. Chimeric
    rows never trigger the error – they are clone-selected and parsed
    per `on_chimeric`. **Exception**: when input is a
    `karyo_preprocessed` object, `"stop"` does not error – unfixable
    rows are returned as NA and a warning is emitted.

  - `"preprocess"`: Apply
    [`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
    to dirty rows. Rows that still have issues after these corrections
    are returned as NA. A message is printed summarising how many rows
    were fixed and how many could not be fixed.

  - `"warn"`: Return NA for all dirty/structural issue rows and emit a
    warning. No fixing is attempted.

  Across all three modes, chimeric rows are clone-selected per
  `on_chimeric` and parsed; a chimeric row becomes NA only when its
  selected clone is itself unparseable (or it has two or more `//`
  separators).

- on_chimeric:

  Character string controlling which clone is parsed for chimeric
  karyotypes (those containing a `//` separator). One of:

  - `"default"`: host clone (before `//`) for normal chimeras; donor
    clone (after `//`) for zero-host chimeras (`.//` prefix).

  - `"host"`: host clone only; zero-host chimeras return NA.

  - `"donor"`: everything after `//` (the donor population). Rows with
    two or more `//` separators (`multiple_chimeric_separator`) are
    always unfixable and return NA. For `karyo_preprocessed` input the
    clone selection was baked in at preprocess time and is inherited;
    passing an explicit, conflicting `on_chimeric` raises an error.

- columns:

  Character vector selecting which output-column classes to include, or
  `NULL` (default) for the full output. `original_karyotype`,
  `preprocessed_karyotype`, `fixable_error`, `unfixable_error`,
  `chimeric_karyotype`, and `chimeric_clone` are always included
  regardless of `columns`. Valid values (any combination):

  - `"classification"`: `normal_karyotype`, `complex_karyotype`,
    `monosomal_karyotype`

  - `"rule"`: rule-matching aberration flags (from `rules`)

  - `"general"`: disease-agnostic structural flags (`general_*`,
    `balanced_translocation`, `unbalanced_translocation`)

  - `"aneuploidy"`: `mono*`/`tris*` monosomy/trisomy flags

  - `"summary"`: `comma_count_aberrations`, `chromosome_count`,
    `total_metaphases`

  - `"der_loss"`: `unbal_partial_loss_<arm>` and `unbal_partial_loss`

- min_metaphases:

  Numeric. Minimum metaphase count for a clone to be eligible to supply
  `chromosome_count`. Default `2`. Clones whose bracket reports fewer
  metaphases are skipped, so a small side clone does not outvote the
  main one; if no clone in a row clears the threshold,
  `chromosome_count` is `NA` for that row. The threshold applies only to
  rows where at least one clone carries a metaphase bracket – a row with
  no brackets anywhere uses all of its clones regardless. Raise it to
  ignore more small clones (`5` is a common cytogenetics convention), or
  set `0` to consider every clone that has a count.

## Value

A tibble with the columns below, grouped by topic rather than listed in
output order:

- original_karyotype: Input karyotype string (always the raw input)

- preprocessed_karyotype: `normalize_iscn()` output of the karyotype
  string, consistent across all `on_issues` modes. In `"preprocess"`
  mode this is also dirty-fixed; in `"warn"`/`"stop"` modes only
  `normalize_iscn()` is applied (chimeric rows are clone-selected in
  every mode). `NA_character_` for issue rows.

- chromosome_count: Integer chromosome count, taken from the most
  abnormal eligible clone in the row (the one whose count is furthest
  from 46). Eligibility is governed by `min_metaphases`; `NA` when no
  clone in the row qualifies, which can happen even though
  `total_metaphases` is non-`NA`.

- One column per aberration flag (0/1 binary); rule-based columns are
  named after a sanitized `flag_name` (see `rules` above)

- monoX, monoY, mono1-22: Monosomy flags for each chromosome

- trisX, trisY, tris1-22: Trisomy flags for each chromosome

- normal_karyotype: 1 if 46,XX or 46,XY, else 0

- total_metaphases: Count from bracket notation

- comma_count_aberrations: Number of comma-separated aberrations

- distinct_aberrations: Number of distinct normalized aberration tokens
  pooled across all clones, excluding sex complements, `idem`/`sl`, and
  range markers. This is the count `complex_karyotype` thresholds at 3.

- complex_karyotype: 1 if \>=3 unique aberrations, else 0

- monosomal_karyotype: 1 if the row has two or more autosomal
  monosomies, or one autosomal monosomy plus at least one structural
  aberration (any `general_*` flag except `general_marker` and
  `general_dmin` – a lone marker or double minute does not count). Two
  carve-outs: sex-chromosome monosomies (`monoX`/`monoY`) never count
  toward either criterion, and a CBF-AML row (`t(8;21)(q22;q22)`,
  `inv(16)(p13q22)`, or `t(16;16)(p13;q22)`) is forced to 0 per clinical
  guidelines even when the criteria are met.

- balanced_translocation: 1 if the row carries a balanced translocation
  – a bare `t(...)` token, or a reciprocal der pair where both partner
  chromosomes appear as centromere donors (e.g.
  `der(5)t(5;17)...,der(17)t(5;17)...`). Independent of
  `general_derivative`. A row may be both balanced and unbalanced.

- unbalanced_translocation: 1 if the row carries an unbalanced
  translocation – a lone der whose reciprocal set is incomplete (e.g.
  `der(a)t(a;b)` or a lone `der(a)t(a;b;c)` three-way der) or a
  whole-arm `der(a;b)`. The specific fusion flag still fires (e.g.
  `der(9)t(9;22)` keeps `t_9_22_q34_q11 = 1`). A der is balanced instead
  only when the full reciprocal set is present (every partner appears as
  a der).

- unbal_partial_loss\_: one column per chromosome arm
  (`unbal_partial_loss_1p`, `unbal_partial_loss_1q`, ...,
  `unbal_partial_loss_22q`, `unbal_partial_loss_Xp`,
  `unbal_partial_loss_Xq`, `unbal_partial_loss_Yp`,
  `unbal_partial_loss_Yq`), set to 1 when an unbalanced der
  translocation implies partial loss of that arm. Kept SEPARATE from
  `del(...)`/`mono*` – an unbalanced-derived 5q loss does not set
  `del_5q`. Derivation needs explicit breakpoints and is limited to
  simple single-junction two-partner `der(a)t(a;b)`; multi-junction
  chains and three-way `t(a;b;c)` derivatives are flagged unbalanced but
  derive no loss. Copy-number-aware gains are out of scope.

- unbal_partial_loss: 1 if an unbalanced der translocation implies any
  partial loss (OR across all `unbal_partial_loss_*` columns).

- fixable_error: 1 if the row had at least one fixable issue (dirty
  marker or chimeric separator) and no unfixable issue, regardless of
  whether the fixable issue was auto-corrected. 0 when
  `unfixable_error = 1`.

- unfixable_error: 1 if row had unfixable structural issues (row is NA),
  else 0

- chimeric_karyotype: 1 if row contained a `//` chimeric separator
  (including zero-host `.//` chimeras and `multiple_chimeric_separator`
  rows), else 0

- chimeric_clone: which clone was parsed for a chimeric row – `"host"`,
  `"donor"`, or `NA_character_` (non-chimeric rows, or chimeric rows
  with no selectable clone such as `multiple_chimeric_separator` or a
  zero-host chimera under `on_chimeric = "host"`)

## See also

[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md),
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md),
[`vignette("data-cleaning")`](https://lucalanino.github.io/karyoparser/articles/data-cleaning.md)

## Examples

``` r
# Basic usage with character vector
karyotypes <- c("46,XX", "47,XY,+21", "46,XX,t(15;17)(q24;q21)")
result <- parse_karyo(karyotypes)

# Usage with data.frame
df <- data.frame(
  sample_id = c("S1", "S2", "S3"),
  iscn = c("46,XX", "47,XY,+21", "46,XX,t(15;17)(q24;q21)")
)
result <- parse_karyo(df, karyotype_column = "iscn")

# Use verbose mode for debugging
result <- parse_karyo(df, karyotype_column = "iscn", verbose = TRUE)
#> Using karyotype column: iscn
#> No ID column given. Use id_column= to carry one through.
#> Parsing 3 karyotype(s).
```
