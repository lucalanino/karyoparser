
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

An R package that turns ISCN karyotype strings into structured binary
features – built for analyzing myeloid neoplasm chromosomal aberrations.

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
  c("46,XX", "47,XY,+21[10]/46,XY[5]", "46,XX,t(9;22)(q34;q11)[20]")
)

result |>
  dplyr::select(original_karyotype, chromosome_count, tris21, t_9_22_q34_q11, normal_karyotype)
#> # A tibble: 3 × 5
#>   original_karyotype     chromosome_count tris21 t_9_22_q34_q11 normal_karyotype
#>   <chr>                             <int>  <int>          <int>            <int>
#> 1 46,XX                                46      0              0                1
#> 2 47,XY,+21[10]/46,XY[5]               47      1              0                0
#> 3 46,XX,t(9;22)(q34;q11…               46      0              1                0
```

## How It Works

| Function / Object | Description |
|----|----|
| `parse_karyo()` | Parse karyotype strings into a wide feature tibble |
| `check_karyo()` | Scan for formatting artifacts and structural errors |
| `preprocess_karyo()` | Clean and normalize dirty strings |
| `validate_rules()` | Build a custom `karyo_rules` object |
| `myeloid_rules` | Default rule set for myeloid neoplasms (data object) |
| `example_karyotypes` | Synthetic karyotypes for trying out the package (data object) |

When working with a new dataset, the recommended pipeline is check,
preprocess, parse: run `check_karyo()` to see what’s wrong,
`preprocess_karyo()` to clean it up, and `parse_karyo()` to extract
features. This way you catch issues before they show up as missing rows.
Passing `on_issues = "preprocess"` to `parse_karyo()` runs the same
check/fix logic internally.

See the [Get
started](https://lucalanino.github.io/karyoparser/articles/karyoparser.html)
article for the full walkthrough: chaining these steps, data frame
input, the `example_karyotypes` dataset, and the complete output column
reference.

## What Is Not Handled

- Sub-band breakpoints: bands are stripped before matching.
- Copy number \> 1: gain/loss is binary; `+8,+8` still just gives
  `tris8 = 1`.
- Partial gain from unbalanced der’s: only the partial loss is flagged
  (`unbal_partial_loss_*`).
- Three-way translocations: only trip the generic
  `general_translocation` flag; partial-loss derivation is limited to
  two-partner `der(a)t(a;b)`.
- Non-ISCN inputs: array CGH/SNP `seq[GRCh38]` notation and standalone
  FISH results aren’t parsed.

## Learn More

- [Get
  started](https://lucalanino.github.io/karyoparser/articles/karyoparser.html)
  – the full pipeline walkthrough, chaining, data frame input, and the
  complete output column reference.
- [Data
  cleaning](https://lucalanino.github.io/karyoparser/articles/data-cleaning.html)
  – the `on_issues` modes, the fixable/unfixable issue catalog, and
  `preprocess_karyo()`’s cleaning pipeline.
- [Chimeric
  karyotypes](https://lucalanino.github.io/karyoparser/articles/chimeric-karyotypes.html)
  – clone selection for `//`-separated chimeric karyotypes.
- [Custom
  rules](https://lucalanino.github.io/karyoparser/articles/custom-rules.html)
  – writing your own rule tables and restricting `parse_karyo()`’s
  output with `columns`.
