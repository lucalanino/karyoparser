# Writing custom rule tables

``` r

library(karyoparser)
```

`myeloid_rules` is
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)’s
default rule set, but it’s just a `karyo_rules` table – a two-column
schema you can build yourself with
[`validate_rules()`](https://lucalanino.github.io/karyoparser/reference/validate_rules.md)
for a different disease area, panel, or one-off lesion.

## The `karyo_rules` schema

A rules table needs exactly two columns:

- `flag_name`: a unique label per rule. Can use ISCN-style punctuation
  for readability (e.g. `"t(9;22)(q34;q11)"`) – the output column is a
  sanitized version, with runs of non-alphanumeric characters replaced
  by `_`.
- `regex`: the pattern matched against each token. Rules fire
  independently – there’s no priority ranking between them, so if two
  patterns can match the same flag, combine them into one `regex` with
  `|` rather than using two rows.

``` r

my_rules <- validate_rules(tibble::tibble(
  flag_name = "t(X;18)(p11;q11)",
  regex = "t\\(X;18\\)\\(p11;q11\\)|t\\(18;X\\)\\(q11;p11\\)"
))
my_rules
#> # A tibble: 1 × 2
#>   flag_name        regex                                              
#> * <chr>            <chr>                                              
#> 1 t(X;18)(p11;q11) "t\\(X;18\\)\\(p11;q11\\)|t\\(18;X\\)\\(q11;p11\\)"

result <- parse_karyo(
  "46,XX,t(X;18)(p11;q11)",
  rules = my_rules,
  on_issues = "warn",
  verbose = FALSE
)
result |>
  dplyr::select(t_X_18_p11_q11, general_translocation)
#> # A tibble: 1 × 2
#>   t_X_18_p11_q11 general_translocation
#>            <int>                 <int>
#> 1              1                     1
```

The output column (`t_X_18_p11_q11`) is free of characters (`(`, `)`,
`;`) that would otherwise force backtick-quoting
(`` r$`t(X;18)(p11;q11)` ``) in interactive use.

## Uniqueness is checked two ways

[`validate_rules()`](https://lucalanino.github.io/karyoparser/reference/validate_rules.md)
rejects duplicate `flag_name` values – each one maps to exactly one
output column:

``` r

validate_rules(tibble::tibble(
  flag_name = c("del(5q)", "del(5q)"),
  regex = c("del\\(5q", "del\\(5\\)\\(q")
))
#> Error:
#> ! `rules` has duplicate `flag_name` value(s): del(5q). Each flag_name must map to exactly one rule -- combine alternative patterns into a single regex instead.
```

It also rejects distinct `flag_name` values that collide once
*sanitized* into the same column name, even though the raw names differ:

``` r

validate_rules(tibble::tibble(
  flag_name = c("t(9;22)", "t(9,22)"),
  regex = c("t\\(9;22\\)", "t\\(9;22\\)")
))
#> Error:
#> ! `rules` has `flag_name` value(s) that collide once sanitized into an output column name: t(9;22), t(9,22). Rename one of the conflicting flag_name values so they remain distinct after stripping non-alphanumeric characters.
```

## Combining rule tables

[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)’s
`rules` argument also accepts a
[`list()`](https://rdrr.io/r/base/list.html) of `karyo_rules` tables,
combining their flags:

``` r

lesion_a <- validate_rules(tibble::tibble(
  flag_name = "flagA",
  regex = "t\\(9;22\\)\\(q34;q11\\)"
))
lesion_b <- validate_rules(tibble::tibble(
  flag_name = "flagB",
  regex = "del\\(5\\)\\(q13q33\\)"
))

result <- parse_karyo(
  c("46,XX,t(9;22)(q34;q11)", "46,XX,del(5)(q13q33)"),
  rules = list(myeloid_rules, lesion_a, lesion_b),
  on_issues = "warn",
  verbose = FALSE
)
result |>
  dplyr::select(t_9_22_q34_q11, flagA, flagB)
#> # A tibble: 2 × 3
#>   t_9_22_q34_q11 flagA flagB
#>            <int> <int> <int>
#> 1              1     1     0
#> 2              0     0     1
```

`flag_name` must stay unique *across* the combined tables too – there’s
no cross-table priority to resolve a clash, so combining two tables that
both define, say, `"shared_flag"` errors instead of silently keeping
one:

``` r

set_a <- validate_rules(tibble::tibble(
  flag_name = "shared_flag",
  regex = "t\\(9;22\\)\\(q34;q11\\)"
))
set_b <- validate_rules(tibble::tibble(
  flag_name = "shared_flag",
  regex = "del\\(5\\)\\(q13q33\\)"
))
parse_karyo(
  "46,XX,t(9;22)(q34;q11)",
  rules = list(set_a, set_b),
  on_issues = "warn",
  verbose = FALSE
)
#> Error:
#> ! `rules` has duplicate `flag_name` value(s): shared_flag. Each flag_name must map to exactly one rule -- combine alternative patterns into a single regex instead.
```

## Rules vs. general flags

General structural categories – translocations, deletions, inversions,
additions, dicentrics, isodicentrics, pseudodicentrics, isochromosomes,
rings, insertions, duplications, triplications, markers, double minutes,
derivatives – are already detected universally by the parser as
`general_*` columns, independently of whatever `rules` table you pass.
You don’t need a rule for “any translocation” – write rules only for the
*specific* lesions your analysis cares about.

## Restricting output with `columns`

The full output is wide – 170 columns for the default `myeloid_rules`.
Pass `columns` to
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
to keep only the groups you need; metadata and status columns
(`original_karyotype`, `fixable_error`, `chimeric_karyotype`, …) are
always included regardless. Valid values: `"classification"`, `"rule"`,
`"general"`, `"aneuploidy"`, `"summary"`, `"der_loss"` (see
[`vignette("karyoparser")`](https://lucalanino.github.io/karyoparser/articles/karyoparser.md)
for what each group contains).

``` r

parse_karyo(
  c("46,XX", "47,XY,+21,t(9;22)(q34;q11)"),
  on_issues = "warn",
  verbose = FALSE,
  columns = c("classification", "general")
)
#> # A tibble: 2 × 26
#>   original_karyotype   preprocessed_karyotype normal_karyotype complex_karyotype
#>   <chr>                <chr>                             <int>             <int>
#> 1 46,XX                46,XX                                 1                 0
#> 2 47,XY,+21,t(9;22)(q… 47,XY,+21,t(9;22)(q34…                0                 0
#> # ℹ 22 more variables: monosomal_karyotype <int>, general_dicentric <int>,
#> #   general_isodicentric <int>, general_isochromosome <int>,
#> #   general_pseudodicentric <int>, general_ring <int>, general_insertion <int>,
#> #   general_duplication <int>, general_triplication <int>,
#> #   general_translocation <int>, general_addition <int>,
#> #   general_inversion <int>, general_deletion <int>, general_marker <int>,
#> #   general_dmin <int>, general_derivative <int>, …
```

`columns = "rule"` reflects whichever rules table you passed – including
a custom one:

``` r

parse_karyo(
  "46,XX,t(X;18)(p11;q11)",
  rules = my_rules,
  on_issues = "warn",
  verbose = FALSE,
  columns = "rule"
)
#> # A tibble: 1 × 7
#>   original_karyotype     preprocessed_karyotype t_X_18_p11_q11 fixable_error
#>   <chr>                  <chr>                           <int>         <int>
#> 1 46,XX,t(X;18)(p11;q11) 46,XX,t(X;18)(p11;q11)              1             0
#> # ℹ 3 more variables: unfixable_error <int>, chimeric_karyotype <int>,
#> #   chimeric_clone <chr>
```
