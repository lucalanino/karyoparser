
<!-- README.md is generated from README.Rmd. Please edit that file -->

# karyoparser <img src="man/figures/logo.png" align="right" height="138" alt="karyoparser logo" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![format-check](https://github.com/lucalanino/karyo-parser/actions/workflows/format-check.yaml/badge.svg)](https://github.com/lucalanino/karyo-parser/actions/workflows/format-check.yaml)
<!-- badges: end -->

An R package for parsing ISCN karyotype strings into structured binary
features. Designed for analysis of myeloid neoplasm-related chromosomal
aberrations.

## Installation

``` r
# install.packages("pak")
pak::pak("lucalanino/karyo-parser")
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

result[, c("original_karyotype", "ploidy_category", "tris21", "t(9;22)(q34;q11)", "normal_karyotype")]
#> # A tibble: 3 × 5
#>   original_karyotype  ploidy_category tris21 `t(9;22)(q34;q11)` normal_karyotype
#>   <chr>               <chr>            <int>              <int>            <int>
#> 1 46,XX               diploid              0                  0                1
#> 2 47,XY,+21[10]/46,X… other                1                  0                0
#> 3 46,XX,t(9;22)(q34;… diploid              0                  1                0
```

One function call, one tibble. Each row is an input karyotype; each
column is a binary flag (0/1), a ploidy category, or a count.

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
result[, c("sample_id", "original_karyotype", "ploidy_category", "tris21")]
#> # A tibble: 3 × 4
#>   sample_id original_karyotype         ploidy_category tris21
#>   <chr>     <chr>                      <chr>            <int>
#> 1 S1        46,XX                      diploid              0
#> 2 S2        47,XY,+21[10]/46,XY[5]     other                1
#> 3 S3        46,XX,t(9;22)(q34;q11)[20] diploid              0
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

| Column | Type | Description |
|----|----|----|
| `original_karyotype` | character | Raw input string |
| `preprocessed_karyotype` | character | Cleaned/normalized string that was parsed; NA for unfixable rows |
| `ploidy_category` | character | `diploid`, `hyperdiploid`, `high_hypodiploid`, `low_hypodiploid`, `near_haploid`, `other`, or `unknown` |
| `chromosome_count` | integer | Count from the most abnormal eligible clone |
| *(aberration flags)* | integer 0/1 | One column per rule – see [Aberration Flags](#aberration-flags) |
| `mono1`–`mono22`, `monoX`, `monoY` | integer 0/1 | Monosomy flags |
| `tris1`–`tris22`, `trisX`, `trisY` | integer 0/1 | Trisomy flags |
| `normal_karyotype` | integer 0/1 | 1 if `46,XX` or `46,XY` exactly |
| `total_metaphases` | integer | Sum of bracket counts; NA if no brackets |
| `comma_count_aberrations` | integer | Aberration count (max across clones; idem-expanded) |
| `complex_karyotype` | integer 0/1 | 1 if \>= 3 unique aberrations |
| `monosomal_karyotype` | integer 0/1 | 1 if \>= 2 autosomal monosomies, or \>= 1 monosomy + \>= 1 structural aberration |
| `mixed_ploidy` | integer 0/1 | 1 if clones span different ploidy categories |
| `fixable_error` | integer 0/1 | 1 if the row had a fixable issue |
| `unfixable_error` | integer 0/1 | 1 if the row had an unfixable issue (row is all NA) |
| `chimeric_karyotype` | integer 0/1 | 1 if input contained a `//` chimeric separator |

## Aberration Flags

All flags output `0` or `1`. Monosomy/trisomy flags are derived from
`-`/`+` tokens; all others come from the rule-matching pipeline.

Rules are organized into **competition groups** (translocation,
inversion, deletion, etc.). Within a group, only the highest-priority
rule fires – so `t(9;22)(q34;q11)` (priority 100) suppresses
`general_translocation` (priority 65). Rules in **different groups
co-fire independently**, so a complex token can contribute to multiple
output flags simultaneously.

| Priority | Rule type |
|----|----|
| 100 | Specific translocations/inversions at canonical breakpoints |
| 95 | Variant-band translocations; isodicentric X |
| 90 | Variable-partner translocations |
| 85 | Chromosome-arm-specific aberrations; isodicentric; pseudodicentric |
| 80 | Dicentric |
| 65–70 | General structural (ring, insertion, duplication, triplication, general translocation) |
| 55–60 | Catch-all general rules (addition, inversion, deletion, marker, derivative) |

**CBF override**: Cases carrying `t(8;21)(q22;q22)`, `inv(16)(p13q22)`,
or `t(16;16)(p13;q22)` are never classified as monosomal, per clinical
guidelines.

Inspect the full default rule set via `myeloid_rules`:

``` r
myeloid_rules
#> # A tibble: 55 × 6
#>    flag_name      regex category priority counts_for_monosomal competition_group
#>  * <chr>          <chr> <chr>       <dbl> <lgl>                <chr>            
#>  1 t(15;17)(q24;… "t\\… specifi…      100 TRUE                 translocation    
#>  2 t(8;21)(q22;q… "t\\… specifi…      100 FALSE                translocation    
#>  3 inv(16)(p13q2… "inv… specifi…      100 FALSE                inversion        
#>  4 t(16;16)(p13;… "t\\… specifi…      100 FALSE                translocation    
#>  5 t(9;11)(p21;q… "t\\… specifi…      100 TRUE                 translocation    
#>  6 t(6;9)(p22;q3… "t\\… specifi…      100 TRUE                 translocation    
#>  7 inv(3)(q21q26) "inv… specifi…      100 TRUE                 inversion        
#>  8 t(3;3)(q21;q2… "t\\… specifi…      100 TRUE                 translocation    
#>  9 t(9;22)(q34;q… "t\\… specifi…      100 TRUE                 translocation    
#> 10 t(1;3)(p36;q2… "t\\… specifi…      100 TRUE                 translocation    
#> # ℹ 45 more rows
```

## Custom Rules

``` r
my_rules <- validate_rules(dplyr::bind_rows(
  as.data.frame(myeloid_rules),
  data.frame(
    flag_name = "t(X;18)(p11;q11)",
    regex = "t\\(X;18\\)\\(p11;q11\\)|t\\(18;X\\)\\(q11;p11\\)",
    category = "specific_tx",
    priority = 100,
    counts_for_monosomal = TRUE,
    competition_group = "translocation"
  )
))

result <- parse_karyo(df, rules = my_rules)
```

Rule columns: `flag_name`, `regex`, `category`, `priority`,
`counts_for_monosomal`, `competition_group`.

## What Is Not Handled

- **Sub-band breakpoints**: Bands are stripped before matching; precise
  breakpoints are not stored.
- **Copy number \> 1**: Gain/loss are binary – `+8,+8` yields
  `tris8 = 1`.
- **Sex chromosome syndromes**: Only simple monosomy/trisomy is captured
  (e.g. no dedicated Turner/Klinefelter flags).
- **Partial monosomy via der/dup**: `der(7)t(1;7)` is not auto-counted
  as monosomy 7.
- **Non-myeloid panels**: Rules are curated for AML/MDS/MPN/CML.
  Lymphoid and solid-tumor lesions are absent.
- **Array CGH / SNP array notation**: `seq[GRCh38]` format is not
  parsed.
- **FISH-only results**: Single-locus FISH outside ISCN strings are not
  processed.
