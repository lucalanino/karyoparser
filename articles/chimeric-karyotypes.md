# Chimeric karyotypes

``` r

library(karyoparser)
```

A `//` in an ISCN karyotype string separates two independent cell
populations – typically host tissue (before `//`) and a donor/engrafted
population (after `//`), as seen post-transplant. karyoparser only
parses one clone per row, so it has to pick a side; the `on_chimeric`
argument (on
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
and
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md))
controls which one.

## `on_chimeric`: choosing a clone

- `"default"`: host clone (before `//`) for normal chimeras; donor clone
  (after `//`) for zero-host chimeras (a leading `.//` prefix – see
  below).
- `"host"`: host clone only; zero-host chimeras return `NA`.
- `"donor"`: everything after `//` (the donor population).

``` r

parse_karyo(
  "47,XY,+8[15]//46,XX[5]",
  on_issues = "warn",
  on_chimeric = "host"
) |>
  dplyr::select(preprocessed_karyotype, chimeric_clone, tris8)
#> Chimeric: 1 (on_chimeric = "host").
#> # A tibble: 1 × 3
#>   preprocessed_karyotype chimeric_clone tris8
#>   <chr>                  <chr>          <int>
#> 1 47,XY,+8[15]           host               1

parse_karyo(
  "46,XX[15]//47,XY,+8[5]",
  on_issues = "warn",
  on_chimeric = "donor"
) |>
  dplyr::select(preprocessed_karyotype, chimeric_clone, tris8)
#> Chimeric: 1 (on_chimeric = "donor").
#> # A tibble: 1 × 3
#>   preprocessed_karyotype chimeric_clone tris8
#>   <chr>                  <chr>          <int>
#> 1 47,XY,+8[5]            donor              1
```

`chimeric_karyotype` is `1` for both rows (any `//` separator), and
`chimeric_clone` records which side was actually parsed.

## Zero-host chimeras

A leading `.//` (or `//`) prefix marks a **zero-host chimera**: a
donor-only population with no host metaphases. Under `"default"` or
`"donor"` the donor clone is parsed as usual; under `"host"` there is no
host clone to parse, so the row is `NA` – but it’s a *policy*-induced
`NA`, not a data-quality problem, so the row still counts as fixable,
not unfixable:

``` r

ck <- check_karyo(".//46,XX,+8[10]")
ck |>
  dplyr::select(karyotype, fixable, unfixable)
#> # A tibble: 1 × 3
#>   karyotype       fixable unfixable
#>   <chr>             <int>     <int>
#> 1 .//46,XX,+8[10]       1         0

parse_karyo(".//46,XX,+8[10]", on_issues = "warn", on_chimeric = "host") |>
  dplyr::select(preprocessed_karyotype, chimeric_clone)
#> Warning: 1 of 1 karyotype(s) had issues and were returned as NA.
#> Chimeric: 1 (on_chimeric = "host").
#> # A tibble: 1 × 2
#>   preprocessed_karyotype chimeric_clone
#>   <chr>                  <chr>         
#> 1 NA                     NA
parse_karyo(".//46,XX,+8[10]", on_issues = "warn", on_chimeric = "default") |>
  dplyr::select(preprocessed_karyotype, chimeric_clone, tris8)
#> Chimeric: 1 (on_chimeric = "default").
#> # A tibble: 1 × 3
#>   preprocessed_karyotype chimeric_clone tris8
#>   <chr>                  <chr>          <int>
#> 1 46,XX,+8[10]           donor              1
```

## Multiple separators are always unfixable

A row with two or more `//` separators (`multiple_chimeric_separator`)
has no well-defined host/donor split, so it’s always unfixable and
returns `NA` regardless of `on_chimeric` or `on_issues`:

``` r

check_karyo("46,XX//47//48") |>
  dplyr::select(karyotype, fixable, unfixable)
#> # A tibble: 1 × 3
#>   karyotype     fixable unfixable
#>   <chr>           <int>     <int>
#> 1 46,XX//47//48       0         1
parse_karyo("46,XX//47//48", on_issues = "preprocess", verbose = FALSE) |>
  dplyr::select(original_karyotype, preprocessed_karyotype, unfixable_error)
#> Chimeric: 1 (on_chimeric = "default").
#> # A tibble: 1 × 3
#>   original_karyotype preprocessed_karyotype unfixable_error
#>   <chr>              <chr>                            <int>
#> 1 46,XX//47//48      NA                                   1
```

## `on_issues` is orthogonal to chimeric handling

Chimeric rows are always clone-selected and parsed, in every `on_issues`
mode – a chimeric row only becomes `NA` when its selected clone is
itself unparseable (or the row has multiple separators). A non-chimeric
issue elsewhere in the batch still triggers `on_issues = "stop"`’s
error; chimeric rows never trigger it themselves:

``` r

parse_karyo(c("46,XX[15]//47,XY,+8[5]", ".46,XX,del(5)(q13)[10]"))
#> Error:
#> ! parse_karyo() stopped: 1 of 2 karyotype(s) have data quality issues.
#>   Fixable:   leading_dot (1)
#> Rerun with on_issues = "preprocess" to auto-correct fixable rows (unfixable rows will be NA).
#> Rerun with on_issues = "warn" to return NA for all issue rows without fixing.
#> Call check_karyo() for a full per-row quality report.
```

``` r

parse_karyo("46,XX[15]//47,XY,+8[5]", verbose = FALSE) |>
  dplyr::select(preprocessed_karyotype, chimeric_clone)
#> Chimeric: 1 (on_chimeric = "default").
#> # A tibble: 1 × 2
#>   preprocessed_karyotype chimeric_clone
#>   <chr>                  <chr>         
#> 1 46,XX[15]              host
```

The second call succeeds under the default `on_issues = "stop"` because
the only issue present is the chimeric separator itself.

## Clone selection and `karyo_preprocessed`

[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
bakes its `on_chimeric` choice into the `preprocessed` column. When that
`karyo_preprocessed` tibble is passed to
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md),
the same choice is inherited automatically – passing an explicit,
conflicting `on_chimeric` at parse time raises an error, since the
donor/host split can’t be un-done after the fact:

``` r

pp <- preprocess_karyo("46,XX[15]//47,XY,+8[5]", on_chimeric = "donor")
parse_karyo(pp) |>
  dplyr::select(preprocessed_karyotype, chimeric_clone, tris8)
#> Chimeric: 1 (on_chimeric = "donor").
#> # A tibble: 1 × 3
#>   preprocessed_karyotype chimeric_clone tris8
#>   <chr>                  <chr>          <int>
#> 1 47,XY,+8[5]            donor              1
```

``` r

parse_karyo(pp, on_chimeric = "host")
#> Error:
#> ! on_chimeric = "host" conflicts with the karyo_preprocessed input, which was created with on_chimeric = "donor". 
#> The donor/host clone selection is baked into the `preprocessed` column, so it cannot be changed here. 
#> Re-run preprocess_karyo(..., on_chimeric = "host"), or pass the raw karyotypes to parse_karyo() with on_chimeric = "host".
```

Re-run `preprocess_karyo(..., on_chimeric = "host")`, or pass the raw
karyotypes straight to `parse_karyo(..., on_chimeric = "host")` instead.
