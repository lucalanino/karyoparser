
<!-- README.md is generated from README.Rmd. Please edit that file -->

# karyoparser <img src="man/figures/logo.png" align="right" height="138" alt="karyoparser logo" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/lucalanino/karyoparser/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/lucalanino/karyoparser/actions/workflows/R-CMD-check.yaml)
[![format-check](https://github.com/lucalanino/karyoparser/actions/workflows/format-check.yaml/badge.svg)](https://github.com/lucalanino/karyoparser/actions/workflows/format-check.yaml)
<!-- badges: end -->

An R package for parsing ISCN karyotype strings into structured binary
features. Designed for analysis of myeloid neoplasm-related chromosomal
aberrations.

## Installation

``` r
# install.packages("pak")
pak::pak("lucalanino/karyoparser")
```

Requires R \>= 4.1.0.

## Quick Start

``` r
library(karyoparser)

result <- parse_karyo(
  c("46,XX", "47,XY,+21[10]/46,XY[5]", "46,XX,t(9;22)(q34;q11)[20]"),
  on_issues = "warn",
  verbose = FALSE
)

result[, c("original_karyotype", "chromosome_count", "tris21", "t(9;22)(q34;q11)", "normal_karyotype")]
#> # A tibble: 3 × 5
#>   original_karyotype chromosome_count tris21 `t(9;22)(q34;q11)` normal_karyotype
#>   <chr>                         <int>  <int>              <int>            <int>
#> 1 46,XX                            46      0                  0                1
#> 2 47,XY,+21[10]/46,…               47      1                  0                0
#> 3 46,XX,t(9;22)(q34…               46      0                  1                0
```

One function call, one tibble. Each row is an input karyotype; each
column is a binary flag (0/1) or a count.

## Main Functions

| Function / Object    | Description                                          |
|----------------------|------------------------------------------------------|
| `parse_karyo()`      | Parse karyotype strings into a wide feature tibble   |
| `check_karyo()`      | Scan for formatting artifacts and structural errors  |
| `preprocess_karyo()` | Clean and normalize dirty strings                    |
| `validate_rules()`   | Build a custom `karyo_rules` object                  |
| `myeloid_rules`      | Default rule set for myeloid neoplasms (data object) |

### Data frame input

`parse_karyo()` accepts a data frame and auto-detects the karyotype
column (from common names: `karyotype`, `iscn`, `karyo`, …) and an
optional id column (`sample_id`, `patient_id`, `id`, `mrn`, …).

``` r
df <- data.frame(
  sample_id = c("S1", "S2", "S3"),
  iscn = c("46,XX", "47,XY,+21[10]/46,XY[5]", "46,XX,t(9;22)(q34;q11)[20]")
)

result <- parse_karyo(df, verbose = FALSE)
result[, c("sample_id", "original_karyotype", "chromosome_count", "tris21")]
#> # A tibble: 3 × 4
#>   sample_id original_karyotype         chromosome_count tris21
#>   <chr>     <chr>                                 <int>  <int>
#> 1 S1        46,XX                                    46      0
#> 2 S2        47,XY,+21[10]/46,XY[5]                   47      1
#> 3 S3        46,XX,t(9;22)(q34;q11)[20]               46      0
```

### `on_issues` parameter

| Value | Dirty rows | Chimeric rows (`//`) | Unfixable rows |
|----|----|----|----|
| `"stop"` (default) | Error immediately | Error immediately | Error immediately |
| `"preprocess"` | Auto-cleaned via `preprocess_karyo()` | Truncated to first clone | Always NA |
| `"warn"` | NA + warning | NA + warning | Always NA |

Use `"stop"` to catch data quality problems early. Switch to
`"preprocess"` once you understand the issues.

## Preprocessing

`preprocess_karyo()` applies a fixed sequence of cleaning steps to dirty
strings:

1.  Trim whitespace
2.  Normalize unicode spaces, dashes, and fullwidth characters
3.  Collapse embedded newlines/tabs
4.  Decode HTML entities (`&lt;`, `&gt;`, `&amp;`)
5.  Strip leading dot(s) before a digit (e.g. `.46,XX` -\> `46,XX`)
6.  Strip FISH/nuc ish suffix
7.  Strip trailing narrative
8.  Collapse mid-string line-wrap artifacts
9.  Insert missing comma after sex chromosome complement
10. Remove space before `mar` token
11. `normalize_iscn()`: whitespace collapsing, delimiter tightening,
    idem/sl/cp normalization

``` r
# Check what issues exist, then clean
ck <- check_karyo(c(".46,XX", "47,XY,+21[10] .Lab note"), verbose = TRUE)
#> Checking 2 karyotype(s)...
#>   Fixable:   2
#>   Unfixable: 0
#>   Fixable breakdown:   leading_dot: 1, trailing_narrative: 1
```

``` r
# Full three-step pipeline -- id column propagates automatically
ck <- check_karyo(df)
clean <- preprocess_karyo(ck)
result <- parse_karyo(clean, verbose = FALSE)
result[, c("sample_id", "original_karyotype", "preprocessed_karyotype")]
#> # A tibble: 3 × 3
#>   sample_id original_karyotype         preprocessed_karyotype    
#>   <chr>     <chr>                      <chr>                     
#> 1 S1        46,XX                      46,XX                     
#> 2 S2        47,XY,+21[10]/46,XY[5]     47,XY,+21[10]/46,XY[5]    
#> 3 S3        46,XX,t(9;22)(q34;q11)[20] 46,XX,t(9;22)(q34;q11)[20]
```

## Output Columns

Each row is one input karyotype; columns fall into a few groups by
provenance. The full output is wide (167 columns for the default
`myeloid_rules`); pass `columns` to `parse_karyo()` to keep only
selected groups – e.g. `columns = c("summary", "balance")` – while
metadata and status columns are always included. Valid values are
`"rule"`, `"general"`, `"aneuploidy"`, `"summary"`, `"balance"`,
`"loss"` (matching the groups below).

| Group | Columns | Type | Description |
|----|----|----|----|
| Metadata | `original_karyotype`, `preprocessed_karyotype` | character | Raw input; cleaned/normalized string that was parsed (NA for unfixable rows) |
| Metadata | `chromosome_count` | integer | Count from the most abnormal eligible clone |
| Rule flags | *(one per `myeloid_rules` entry)* | integer 0/1 | Specific lesions – see [Aberration Flags](#aberration-flags) |
| General flags | `general_translocation`, `general_deletion`, `general_inversion`, `general_addition`, `general_dicentric`, `general_isodicentric`, `general_pseudodicentric`, `general_isochromosome`, `general_ring`, `general_insertion`, `general_duplication`, `general_triplication`, `general_marker`, `general_derivative` | integer 0/1 | Universal structural-aberration detections |
| Aneuploidy | `mono1`–`mono22`, `monoX`, `monoY`; `tris1`–`tris22`, `trisX`, `trisY` | integer 0/1 | Whole-chromosome loss/gain from `-`/`+` tokens |
| Summary | `normal_karyotype` | integer 0/1 | 1 if `46,XX` or `46,XY` exactly |
| Summary | `comma_count_aberrations` | integer | Aberration count (max across clones; idem-expanded) |
| Summary | `complex_karyotype` | integer 0/1 | 1 if \>= 3 distinct aberrations across clones |
| Summary | `monosomal_karyotype` | integer 0/1 | 1 if \>= 2 autosomal monosomies, or \>= 1 monosomy + \>= 1 structural aberration |
| Translocation balance | `balanced_translocation`, `unbalanced_translocation` | integer 0/1 | Whether the row carries balanced / unbalanced translocations |
| Derived loss | `unbal_partial_loss_<arm>` (one per chromosome arm), `unbal_partial_loss` | integer 0/1 | Partial arm loss implied by an unbalanced der translocation |
| Status | `total_metaphases` | integer | Sum of bracket counts; NA if no brackets |
| Status | `fixable_error`, `unfixable_error` | integer 0/1 | Row had a fixable / unfixable issue (unfixable rows are all NA) |
| Status | `chimeric_karyotype` | integer 0/1 | 1 if input contained a `//` chimeric separator |
| Status | `chimeric_clone` | character | Which clone was parsed for a chimeric row (`"host"`, `"donor"`, or NA) |

### Return object classes and attributes

`check_karyo()` and `preprocess_karyo()` return tibbles with an extra S3
class, so their output can be piped straight into the next pipeline step
without re-specifying arguments:

| Function | Class | Notes |
|----|----|----|
| `check_karyo()` | `karyo_check` (+ `tbl_df`) | Accepted as `preprocess_karyo()` input; the assessment is reused instead of re-run. |
| `preprocess_karyo()` | `karyo_preprocessed` (+ `tbl_df`) | Accepted as `parse_karyo()` input; the `preprocessed` column and cached id/issue info are reused automatically. |
| `parse_karyo()` | plain `tbl_df` | No custom class. |

`parse_karyo()` output also carries a `karyoparser_version` attribute
(`attr(result, "karyoparser_version")`), set to the installed package
version, so a saved result can be traced back to the version that produced
it.

## Aberration Flags

All flags output `0` or `1`, and **every flag fires independently** –
there is no priority or competition between them, so a single token can
light up several flags at once. Flags come from three sources:

- **Specific rule flags** – one per `myeloid_rules` entry, fired when
  its regex matches a token (e.g. `t(9;22)(q34;q11)`, `del(5q)`).
- **General structural flags** (`general_*`) – universal,
  disease-agnostic detections (any translocation, deletion, inversion,
  addition, dicentric, isodicentric, pseudodicentric, isochromosome,
  ring, insertion, duplication, triplication, marker, derivative).
  Always computed, independent of the rule set.
- **Aneuploidy flags** (`mono*`/`tris*`) – whole-chromosome loss/gain
  from `-`/`+` tokens.

Because flags are independent, a token contributes to **all** matching
flags. For example `t(9;11)(p21;q23)` sets the specific
`t(9;11)(p21;q23)`, the family flag `t(v;11q23)`, and
`general_translocation`; `del(5q)` sets both `del(5q)` and
`general_deletion`.

- **Family detectors** (`t(v;11q23)`, `t(v;11p15)`, `t(3q26;v)`,
  `t(5q)`, `t(12p)`) deliberately co-fire with the more specific flags –
  they mean “any rearrangement involving this region/arm”.
- **`*_other` variants** (e.g. `t(9;22)_other`) are the complement of
  the canonical breakpoints: they fire for the same chromosome pair at
  *non-canonical* breakpoints only (the exclusion is encoded in the
  regex, so a canonical token sets only the specific flag).

**CBF override**: Cases carrying `t(8;21)(q22;q22)`, `inv(16)(p13q22)`,
or `t(16;16)(p13;q22)` are never classified as monosomal, per clinical
guidelines.

Inspect the full default rule set via `myeloid_rules`:

``` r
myeloid_rules
#> # A tibble: 42 × 3
#>    flag_name         regex                                              category
#>  * <chr>             <chr>                                              <chr>   
#>  1 t(15;17)(q24;q21) "t\\(15;17\\)\\((q24|q22);q21\\)|t\\(17;15\\)\\(q… specifi…
#>  2 t(8;21)(q22;q22)  "t\\(8;21\\)\\((q22|q21(\\.3)?);q22\\)|t\\(21;8\\… specifi…
#>  3 inv(16)(p13q22)   "inv\\(16\\)\\(p13q22\\)"                          specifi…
#>  4 t(16;16)(p13;q22) "t\\(16;16\\)\\(p13;q22\\)"                        specifi…
#>  5 t(9;11)(p21;q23)  "t\\(9;11\\)\\(p21;q23\\)|t\\(11;9\\)\\(q23;p21\\… specifi…
#>  6 t(6;9)(p22;q34)   "t\\(6;9\\)\\(p22;q34\\)|t\\(9;6\\)\\(q34;p22\\)"  specifi…
#>  7 inv(3)(q21q26)    "inv\\(3\\)\\(q21q26\\)"                           specifi…
#>  8 t(3;3)(q21;q26)   "t\\(3;3\\)\\(q21;q26\\)"                          specifi…
#>  9 t(9;22)(q34;q11)  "t\\(9;22\\)\\(q34;q11\\)|t\\(22;9\\)\\(q11;q34\\… specifi…
#> 10 t(1;3)(p36;q21)   "t\\(1;3\\)\\(p36;q21\\)|t\\(3;1\\)\\(q21;p36\\)"  specifi…
#> # ℹ 32 more rows
```

## Custom Rules

``` r
my_rules <- validate_rules(dplyr::bind_rows(
  as.data.frame(myeloid_rules),
  data.frame(
    flag_name = "t(X;18)(p11;q11)",
    regex = "t\\(X;18\\)\\(p11;q11\\)|t\\(18;X\\)\\(q11;p11\\)",
    category = "specific_tx"
  )
))

result <- parse_karyo(df, rules = my_rules)
```

Rule columns: `flag_name`, `regex`, `category`. Each rule fires
independently when its `regex` matches a token – there is no priority
ranking, so if you need a flag to *not* fire in some case (e.g. a
“non-canonical breakpoints only” variant), encode that exclusion in the
`regex` (e.g. with a negative lookahead). General structural categories
(translocations, deletions, dicentrics, …) are detected universally by
the parser and do **not** need to be added as rules.

## What Is Not Handled

- **Sub-band breakpoints**: Bands are stripped before matching; precise
  breakpoints are not stored.
- **Copy number \> 1**: Gain/loss are binary – `+8,+8` yields
  `tris8 = 1`.
- **Partial gain, and offsetting of partial loss**: An unbalanced
  `der(a)t(a;b)` records the implied partial *loss*
  (`unbal_partial_loss_*`) but never a partial *gain* – there is no
  `unbal_partial_gain_*`. The loss is derived token-locally from the der
  breakpoints: it fires purely from the der’s structure and is **not**
  reconciled against the rest of the karyotype. So a co-occurring gain
  of the same material (e.g. `+9` or a `dup`) that would offset the loss
  does **not** suppress it – `47,XX,+9,der(9)t(9;22)(q34;q11)` still
  reports `unbal_partial_loss_9q = 1` even though distal 9q sits at two
  copies.
- **Sex chromosome syndromes**: Only simple monosomy/trisomy is captured
  (e.g. no dedicated Turner/Klinefelter flags).
- **Partial monosomy via der/dup**: `der(7)t(1;7)` is not auto-counted
  as monosomy 7.
- **`idem` token-level expansion**: `idem` (stemline shorthand) is
  expanded for aberration *counting* but not for token-level
  flag/balance analysis. Consequence: a reciprocal der pair split across
  the `idem` boundary (stemline carries `der(a)`, an `idem` subclone
  adds `der(b)`) is reported `unbalanced_translocation = 1` only,
  missing the `balanced_translocation = 1` the idem-inherited partner
  implies.
- **Multi-partner (three-way) translocations**: A complete `t(a;b;c)`
  matches only the generic `general_translocation` flag (no
  specific-breakpoint flag) and is assumed balanced. A `der()` of a
  three-way is classified for translocation balance – a lone der is
  `unbalanced_translocation`, a complete der-by-der reciprocal set is
  `balanced_translocation` – but derives no partial loss
  (`unbal_partial_loss_*`); partial-loss derivation remains limited to
  simple two-partner `der(a)t(a;b)`. Both still count toward
  `complex_karyotype`/`monosomal_karyotype`.
- **Non-myeloid panels**: Rules are curated for AML/MDS/MPN/CML.
  Lymphoid and solid-tumor lesions are absent.
- **Mosaicism (`mos`)**: Karyotypes with a leading `mos` prefix
  (e.g. `mos 47,XXY[10]/46,XY[5]`) are flagged `mosaic_karyotype` and
  left unparsed.
- **Non-clonal single-cell abnormalities (`ncSCA`)**: Karyotypes
  carrying an `ncSCA` token are flagged `non_clonal_sca` and left
  unparsed.
- **Constitutional abnormalities (`c`)**: Sex complements with a
  constitutional `c` suffix (e.g. `47,XXYc`) are flagged
  `constitutional_sex_complement` rather than parsed as acquired.
- **Array CGH / SNP array notation**: `seq[GRCh38]` format is not
  parsed.
- **FISH-only results**: Single-locus FISH outside ISCN strings are not
  processed.
