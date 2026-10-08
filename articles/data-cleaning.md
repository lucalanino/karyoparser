# Cleaning dirty karyotype strings

``` r

library(karyoparser)
```

Real-world karyotype strings are messy: stray whitespace, HTML entities,
missing commas, FISH annotations tacked onto the end, and outright
out-of-scope constructs like mosaicism. This vignette covers
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md),
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md),
and the `on_issues` argument to
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
in detail.

## Fixable vs. unfixable

Every issue
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md)
finds is either **fixable** (a formatting artifact) or **unfixable** (a
structural problem).
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
only touches the fixable kind – unfixable rows always come out `NA`, no
matter which `on_issues` mode you use.

``` r

ck <- check_karyo(example_karyotypes, karyotype_column = "karyotype")
names(ck)
#>  [1] "karyotype"                     "fixable"                      
#>  [3] "unfixable"                     "unfixable_reason"             
#>  [5] "chimeric_no_count"             "constitutional_sex_complement"
#>  [7] "empty"                         "invalid_breakpoint"           
#>  [9] "invalid_idem"                  "mosaic_karyotype"             
#> [11] "multiple_chimeric_separator"   "no_chromosome_count"          
#> [13] "no_sex_complement"             "single_token"                 
#> [15] "stray_non_ascii"               "unbalanced_brackets"          
#> [17] "unbalanced_parentheses"        "unparseable_bracket"          
#> [19] "unrecognized_gain_loss"        "unrecognized_token"           
#> [21] "updated_iscn"
```

`unfixable_reason` spells out why a row is unfixable, naming the
offending token where there is one
(e.g. `unrecognized_token: Cannot classify 'idem*2' as an ISCN token`).
The columns after it are one 0/1 flag per *unfixable* issue type
(alphabetical). Fixable issue types don’t get their own columns –
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
resolves them automatically, so a row is simply flagged `fixable`,
without needing to know which specific fixable pattern fired. The full
catalog of issue types, whether or not they surface as a column:

- **Fixable**: `breakpoint_semicolon` (`t(9;16)(q34p13)` for
  `t(9;16)(q34;p13)`), `case_notation` (a lowercase sex chromosome such
  as `46,xy`, `-y` or `t(x;5)`, or an uppercase ISCN keyword such as
  `Del(5)`, `MAR` or `(Q13Q33)`; a lowercase `x` copy multiplier like
  `der(1)x2` is left alone), `chimeric_separator`, `colon_separator`
  (`t(5:17)(q13:q11)`), `count_sex_separator` (e.g. `46XY` or `45.XY`
  missing the comma before the sex complement), `dot_separator`
  (`del(5)(q13q33).-7`), `embedded_newline`, `fish_notation`,
  `fullwidth_punctuation` (fullwidth brackets, tilde, parentheses,
  comma, semicolon, or equals sign standing in for their ASCII
  equivalents), `gain_loss_comma` (`-7+mar`), `html_entities`,
  `iso_indicator` (non-ISCN `iso(17q)` written for `i(17q)`),
  `leading_dot`, `mar_space`, `midstring_linewrap`, `missing_sex_comma`
  (e.g. `46,XX der(...)` or `46,XY+8` missing the comma after the sex
  complement – a letter-glued form like `46,XYdel(5q)` is auto-fixed
  when the glued token is a recognized aberration indicator, otherwise
  it’s detected but comes back `NA`), `non_ascii_homoglyph` (a Greek
  Chi/Upsilon standing in for Latin `X`/`Y` in the sex complement),
  `non_clonal_sca` (an `ncSCA` token, standalone or parenthesized),
  `paren_comma` (`del(5)(q11,q33)`; a decimal comma `q11,2` becomes
  `q11.2`), `star_multiplier` (`+mar*2` for `+marx2`),
  `trailing_narrative`, `underscore_prefix` (`46,_XY,_del(5)...`),
  `unicode_notation`, `zero_host_chimera`.
- **Unfixable**: `chimeric_no_count`, `constitutional_sex_complement`,
  `empty`, `invalid_breakpoint`, `invalid_idem`, `mosaic_karyotype`,
  `multiple_chimeric_separator`, `no_chromosome_count`,
  `no_sex_complement`, `single_token`, `stray_non_ascii`,
  `unbalanced_brackets`, `unbalanced_parentheses`,
  `unparseable_bracket`, `unrecognized_gain_loss`, `unrecognized_token`,
  `updated_iscn`.

## Unfixable issues: what to do

Unfixable rows always come back `NA`, but the *reason* varies – some are
genuine data problems worth chasing down in the source record, others
are constructs the package deliberately doesn’t handle.

| Issue type | What it usually means | What to do |
|----|----|----|
| `empty` | Input was `NA` or `""` | Confirm the source column doesn’t have missing/blank karyotype values before parsing, or filter them out upstream. |
| `no_chromosome_count` | String, or any later clone, doesn’t start with a chromosome count (e.g. `46,XY[20]/text:46,idem,...`) | Check whether a leading chromosome count was stripped upstream (e.g. spreadsheet auto-formatting), or whether the string is karyotype notation at all. |
| `no_sex_complement` | First clone has no sex chromosome complement (`XX`/`XY`/etc.) | Check whether the string was truncated before the complement, uses non-standard complement notation, or (if `missing_sex_comma` is also flagged) has the complement glued directly to an unrecognized or garbled aberration token with no comma, e.g. `47,XY(inv)(9)`. A glued *recognized* aberration indicator, e.g. `46,XYdel(5q)`, is auto-fixed and won’t reach this check. |
| `single_token` | Bare comma-less string, e.g. `"8"` | Ambiguous between a chromosome count and an abnormality whose leading `+`/`-` sign was stripped – a common artifact of opening ISCN strings in Excel. Check the original cell format/source. |
| `invalid_idem` | `idem` appears in the first clone | `idem` means “same as the previous clone”, but there is no previous clone here. Check whether an earlier clone was dropped from the string. |
| `unbalanced_parentheses` / `unbalanced_brackets` | Mismatched `(`/`)` or `[`/`]` counts | Usually a truncation or copy-paste error. Inspect the original source record. |
| `unparseable_bracket` | Bracket content isn’t a valid metaphase count (`[n]`, `[cpN]`, `[n~m]`) | Look for stray text inside the brackets, e.g. a note or annotation that ended up bracketed. |
| `unrecognized_gain_loss` | A `+`/`-` token that is neither a whole-chromosome gain/loss (`+8`, `-Xx2`), an unidentified ring/marker/double minute (`+1~3r`, `+2mar`), nor a structural aberration (`+der(1)...`), e.g. `+8q` or `+23` | Monosomy/trisomy flags only fire on whole tokens, so an unforeseen shape would otherwise be silently skipped. Correct the token in the source record. |
| `unrecognized_token` | Any other aberration token with no recognized ISCN shape, e.g. free text, a bare number (`46,XY,8,del(5)(q13q33)[20]`), or an indicator with no breakpoints (`47,XY,+8,del[20]`) | Unrecognized tokens would otherwise be counted as aberrations and inflate `complex_karyotype`. Only token structure is checked (indicator, parenthesized groups, suffixes), never band contents. Correct the token in the source record. |
| `invalid_breakpoint` | A recognized token whose parentheses hold an impossible chromosome (`der(+)`, `inv(1198q15q22)`) or breakpoint: an empty entry (`()`, `(p35;;p24)`), a band without an arm (`(q12;21)`), or a breakpoint count that doesn’t match the chromosomes (`t(5;2)(q11)`) | ISCN uncertainty (`?`, `q1?`), `pter`/`qter`/`cen`, `or` alternatives and uncertain ranges (`q13-14`, `q33~34`) are accepted. Correct the breakpoints in the source record. |
| `multiple_chimeric_separator` | Two or more `//` separators | The parser can’t auto-resolve more than one chimeric boundary. Decide manually which clone population to keep. |
| `chimeric_no_count` | A `//` not followed by a chromosome count, e.g. `46,XX[10]//XY[5]` | Every population after `//` needs its own ploidy. A missing count is never guessed to be a `/` typo; check the source record. |
| `updated_iscn` | Contains an `"Updated ISCN"` marker | The string may have been superseded by a later correction. Check the source record for the corrected karyotype. |
| `constitutional_sex_complement` | Sex complement has a `c` suffix, e.g. `47,XXYc` | Out of scope by design, not a data-quality issue – see [Out-of-scope constructs](#out-of-scope-constructs-take-precedence) below. |
| `mosaic_karyotype` | Leading `mos` prefix | Out of scope by design, not a data-quality issue – see [Out-of-scope constructs](#out-of-scope-constructs-take-precedence) below. |
| `stray_non_ascii` | Leftover non-ASCII character with no known automatic fix | Identify and manually normalize or remove the character. If it’s a recurring pattern in your data, consider filing an issue so it can be added as an automatic fix. |

## Out-of-scope constructs take precedence

Mosaicism (`mos`) and constitutional abnormalities (a sex complement
with a `c` suffix, e.g. `47,XXYc`) are deliberately out of scope and
flagged unfixable – even though the string may otherwise look
structurally complete. These take precedence over
`no_chromosome_count`/`no_sex_complement`, so a present count or
complement isn’t mislabeled as missing:

``` r

ck <- check_karyo(c(
  "mos 47,XXY[10]/46,XY[5]",
  "47,XXYc[20]"
))
ck |>
  dplyr::select(
    karyotype,
    mosaic_karyotype,
    constitutional_sex_complement,
    no_chromosome_count,
    no_sex_complement
  )
#> # A tibble: 2 × 5
#>   karyotype          mosaic_karyotype constitutional_sex_c…¹ no_chromosome_count
#>   <chr>                         <int>                  <int>               <int>
#> 1 mos 47,XXY[10]/46…                1                      0                   0
#> 2 47,XXYc[20]                       0                      1                   0
#> # ℹ abbreviated name: ¹​constitutional_sex_complement
#> # ℹ 1 more variable: no_sex_complement <int>
```

Each row has a chromosome count and sex complement present, but
`no_chromosome_count`/`no_sex_complement` stay `0` – the out-of-scope
flag fires instead.

## `ncSCA` tokens are silently stripped

Non-clonal single-cell abnormalities (an `ncSCA` token, standalone or
parenthesized) carry no structural detail of their own, so
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
strips the token and parses whatever clone remains:

``` r

preprocess_karyo(c("ncSCA[4]/46,XY[11]", "46,XX(ncSCA)[1]//46,XY[19]"))
#> # A tibble: 2 × 3
#>   original                   preprocessed status
#>   <chr>                      <chr>        <chr> 
#> 1 ncSCA[4]/46,XY[11]         46,XY[11]    fixed 
#> 2 46,XX(ncSCA)[1]//46,XY[19] 46,XX[1]     fixed
```

## `preprocess_karyo()`’s cleaning pipeline

[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
is the only function that rewrites a string. It applies fixes in a fixed
order, then a final `normalize_iscn()` pass:

1.  Trim leading/trailing whitespace
2.  Replace unicode spaces (NBSP), zero-width spaces, dashes
    (em/en-dash), fullwidth `+`/`-` characters
3.  Normalize Greek letter homoglyphs (Chi/Upsilon) to Latin `X`/`Y` in
    the sex chromosome complement
4.  Replace fullwidth ISCN punctuation (brackets, tilde, parentheses,
    comma, semicolon, equals sign) with their ASCII equivalents
5.  Collapse embedded newlines/tabs to spaces
6.  Decode HTML entities (`&lt;` -\> `<`, `&gt;` -\> `>`, `&amp;` -\>
    `&`)
7.  Strip leading dot(s) before a digit (e.g. `.46,XX` -\> `46,XX`)
8.  Strip underscore token prefixes (`46,_XY,_-7` -\> `46,XY,-7`)
9.  Canonicalize letter case: uppercase sex chromosomes (`46,xy` -\>
    `46,XY`, `-y` -\> `-Y`), lowercase ISCN keywords (`Del(5)(Q13Q33)`
    -\> `del(5)(q13q33)`)
10. Rewrite the non-ISCN `iso(` isochromosome indicator to `i(`
11. Insert missing/dotted separator between chromosome count and sex
    complement (e.g. `46XY` or `45.XY` -\> `46,XY`/`45,XY`)
12. Strip FISH/`nuc ish` annotation (suffix or mid-clone before
    metaphase count)
13. Strip trailing narrative: `] .text` -\> `]`; `) Capital text` -\>
    `)`; a trailing `<AML>`-style annotation
14. Collapse mid-string line-wrap artifacts (`, .der(...)` or
    `,.der(...)` -\> `,der(...)`)
15. Insert missing comma after sex chromosome complement
16. Remove space between count and `mar` token (`+1~4 mar` -\>
    `+1~4mar`)
17. Strip non-clonal single-cell abnormality (`ncSCA`) tokens,
    standalone or parenthesized
18. Repair separators: a comma inside parentheses -\> `;`
    (`del(5)(q11,q33)`); a colon used for `;` (`t(5:17)(q13:q11)`),
    leaving the detailed ISCN system (`::`, `->`) alone; a missing `;`
    between two chromosomes’ breakpoints (`t(9;16)(q34p13)` -\>
    `t(9;16)(q34;p13)`); a dot used for `,` between tokens
    (`del(5)(q13q33).-7`, `45.idem`); `*N` -\> `xN`; a missing comma
    before a `+`/`-` token (`-7+mar` -\> `-7,+mar`)
19. `normalize_iscn()`: whitespace collapsing, delimiter tightening,
    idem/sl/cp normalization

(`zero_host_chimera` and chimeric clone selection happen alongside this
pipeline but are covered separately in
[`vignette("chimeric-karyotypes")`](https://lucalanino.github.io/karyoparser/articles/chimeric-karyotypes.md).)

``` r

messy <- c(
  ".46,XX,del(5)(q13)[10]",
  "46,XX,t(9;22)(q34;q11)[15] &lt;AML&gt;",
  "46,XX[cp 20]"
)
preprocess_karyo(messy, verbose = TRUE)
#> Preprocessing 3 karyotype(s)...
#>   Fixed:      2
#>   Clean:      1
#> # A tibble: 3 × 3
#>   original                               preprocessed               status
#>   <chr>                                  <chr>                      <chr> 
#> 1 .46,XX,del(5)(q13)[10]                 46,XX,del(5)(q13)[10]      fixed 
#> 2 46,XX,t(9;22)(q34;q11)[15] &lt;AML&gt; 46,XX,t(9;22)(q34;q11)[15] fixed 
#> 3 46,XX[cp 20]                           46,XX[cp20]                clean
```

[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
accepts a `karyo_check` object directly (as returned by
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md)),
reusing its cached assessment instead of re-scanning:

``` r

ck <- check_karyo(example_karyotypes, karyotype_column = "karyotype")
pp <- preprocess_karyo(ck)
pp
#> # A tibble: 102 × 3
#>    original                    preprocessed                status
#>    <chr>                       <chr>                       <chr> 
#>  1 46,XX[20]                   46,XX[20]                   clean 
#>  2 46,XY[20]                   46,XY[20]                   clean 
#>  3 46,XY                       46,XY                       clean 
#>  4 46,XX,t(15;17)(q24;q21)[20] 46,XX,t(15;17)(q24;q21)[20] clean 
#>  5 46,XY,t(8;21)(q22;q22)[18]  46,XY,t(8;21)(q22;q22)[18]  clean 
#>  6 46,XX,inv(16)(p13q22)[15]   46,XX,inv(16)(p13q22)[15]   clean 
#>  7 46,XY,t(16;16)(p13;q22)[12] 46,XY,t(16;16)(p13;q22)[12] clean 
#>  8 46,XX,t(9;11)(p21;q23)[14]  46,XX,t(9;11)(p21;q23)[14]  clean 
#>  9 46,XY,t(6;9)(p22;q34)[16]   46,XY,t(6;9)(p22;q34)[16]   clean 
#> 10 46,XX,inv(3)(q21q26)[10]    46,XX,inv(3)(q21q26)[10]    clean 
#> # ℹ 92 more rows
```

## `on_issues`: how `parse_karyo()` handles issues

[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)’s
`on_issues` argument controls what happens to dirty and structural rows.
It has no effect on chimeric handling – see
[`vignette("chimeric-karyotypes")`](https://lucalanino.github.io/karyoparser/articles/chimeric-karyotypes.md)
– chimeric rows are always clone-selected and parsed regardless of
`on_issues`; only a non-chimeric issue on the *same* row prevents that.

### `"stop"` (default)

Raises an error immediately if any non-chimeric issue is found, listing
fixable/unfixable issue types and counts:

``` r

parse_karyo(c("46,XX", ".46,XX,del(5)(q13)[10]", "mos 47,XXY[10]"))
#> Error:
#> ! parse_karyo() stopped: 2 of 3 karyotype(s) have data quality issues.
#>   Fixable:   leading_dot (1)
#>   Unfixable: mosaic_karyotype (1)
#> Rerun with on_issues = "preprocess" to auto-correct fixable rows (unfixable rows will be NA).
#> Rerun with on_issues = "warn" to return NA for all issue rows without fixing.
#> Call check_karyo() for a full per-row quality report.
```

A non-chimeric issue anywhere in the batch triggers this error even if
other rows in the same batch are chimeric – chimeric rows themselves
never trigger it.

### `"preprocess"`

Applies
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
to dirty rows; rows still broken afterwards (the unfixable kind) come
back `NA`:

``` r

parse_karyo(
  c("46,XX", ".46,XX,del(5)(q13)[10]", "mos 47,XXY[10]"),
  on_issues = "preprocess",
  verbose = TRUE
)
#>   Fixed:      1
#>   Unfixable:  1  (returned as NA)
#> Parsing 2 karyotype(s).
#> # A tibble: 3 × 170
#>   original_karyotype   preprocessed_karyotype normal_karyotype complex_karyotype
#>   <chr>                <chr>                             <int>             <int>
#> 1 46,XX                46,XX                                 1                 0
#> 2 .46,XX,del(5)(q13)[… 46,XX,del(5)(q13)[10]                 0                 0
#> 3 mos 47,XXY[10]       NA                                   NA                NA
#> # ℹ 166 more variables: monosomal_karyotype <int>, t_15_17_q24_q21 <int>,
#> #   t_8_21_q22_q22 <int>, inv_16_p13q22 <int>, t_16_16_p13_q22 <int>,
#> #   t_9_11_p21_q23 <int>, t_6_9_p22_q34 <int>, inv_3_q21q26 <int>,
#> #   t_3_3_q21_q26 <int>, t_9_22_q34_q11 <int>, t_1_3_p36_q21 <int>,
#> #   t_1_22_p13_q13 <int>, t_3_5_q25_q35 <int>, t_5_11_q35_p15 <int>,
#> #   t_7_12_q36_p13 <int>, t_8_16_p11_p13 <int>, t_10_11_p12_q14 <int>,
#> #   t_11_12_p15_p13 <int>, t_16_21_p11_q22 <int>, t_16_21_q24_q22 <int>, …
```

### `"warn"`

Returns `NA` for all dirty/structural issue rows and emits a warning –
no fixing is attempted:

``` r

parse_karyo(
  c("46,XX", ".46,XX,del(5)(q13)[10]", "mos 47,XXY[10]"),
  on_issues = "warn"
)
#> Warning: 2 of 3 karyotype(s) had issues and were returned as NA.
#> # A tibble: 3 × 170
#>   original_karyotype   preprocessed_karyotype normal_karyotype complex_karyotype
#>   <chr>                <chr>                             <int>             <int>
#> 1 46,XX                46,XX                                 1                 0
#> 2 .46,XX,del(5)(q13)[… NA                                   NA                NA
#> 3 mos 47,XXY[10]       NA                                   NA                NA
#> # ℹ 166 more variables: monosomal_karyotype <int>, t_15_17_q24_q21 <int>,
#> #   t_8_21_q22_q22 <int>, inv_16_p13q22 <int>, t_16_16_p13_q22 <int>,
#> #   t_9_11_p21_q23 <int>, t_6_9_p22_q34 <int>, inv_3_q21q26 <int>,
#> #   t_3_3_q21_q26 <int>, t_9_22_q34_q11 <int>, t_1_3_p36_q21 <int>,
#> #   t_1_22_p13_q13 <int>, t_3_5_q25_q35 <int>, t_5_11_q35_p15 <int>,
#> #   t_7_12_q36_p13 <int>, t_8_16_p11_p13 <int>, t_10_11_p12_q14 <int>,
#> #   t_11_12_p15_p13 <int>, t_16_21_p11_q22 <int>, t_16_21_q24_q22 <int>, …
```

## Next steps

- [`vignette("chimeric-karyotypes")`](https://lucalanino.github.io/karyoparser/articles/chimeric-karyotypes.md)
  – clone selection for `//`-separated chimeric karyotypes.
- [`vignette("custom-rules")`](https://lucalanino.github.io/karyoparser/articles/custom-rules.md)
  – writing your own rule tables.
