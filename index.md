# karyoparser

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
| [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md) | Parse karyotype strings into a wide feature tibble |
| [`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md) | Scan for formatting artifacts and structural errors |
| [`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md) | Clean and normalize dirty strings |
| [`assign_risk()`](https://lucalanino.github.io/karyoparser/reference/assign_risk.md) | Assign cytogenetic risk categories to a [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md) result |
| [`validate_rules()`](https://lucalanino.github.io/karyoparser/reference/validate_rules.md) | Build a custom `karyo_rules` object |
| `myeloid_rules` | Default rule set for myeloid neoplasms (data object) |
| `example_karyotypes` | Synthetic karyotypes for trying out the package (data object) |

When working with a new dataset, the recommended pipeline is check,
preprocess, parse: run
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md)
to see what’s wrong,
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
to clean it up, and
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
to extract features. This way you catch issues before they show up as
missing rows. Passing `on_issues = "preprocess"` to
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
runs the same check/fix logic internally.

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
  [`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)’s
  cleaning pipeline.
- [Chimeric
  karyotypes](https://lucalanino.github.io/karyoparser/articles/chimeric-karyotypes.html)
  – clone selection for `//`-separated chimeric karyotypes.
- [Risk
  stratification](https://lucalanino.github.io/karyoparser/articles/risk-stratification.html)
  – IPSS-R and ELN 2022 cytogenetic risk categories via
  [`assign_risk()`](https://lucalanino.github.io/karyoparser/reference/assign_risk.md).
- [Custom
  rules](https://lucalanino.github.io/karyoparser/articles/custom-rules.html)
  – writing your own rule tables and restricting
  [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)’s
  output with `columns`.
