# Getting started with karyoparser

``` r

library(karyoparser)
```

karyoparser turns ISCN karyotype strings into a wide tibble of binary
features. This vignette walks through the three pipeline functions, the
`example_karyotypes` dataset that ships with the package, and how to
read the output.

## The three pipeline functions

| Function / Object | Description |
|----|----|
| [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md) | Parse karyotype strings into a wide feature tibble |
| [`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md) | Scan for formatting artifacts and structural errors |
| [`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md) | Clean and normalize dirty strings |
| [`assign_risk()`](https://lucalanino.github.io/karyoparser/reference/assign_risk.md) | Assign cytogenetic risk categories to a [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md) result |
| [`validate_rules()`](https://lucalanino.github.io/karyoparser/reference/validate_rules.md) | Build a custom `karyo_rules` object |
| `myeloid_rules` | Default rule set for myeloid neoplasms (data object) |
| `example_karyotypes` | Synthetic karyotypes for trying out the package (data object) |

[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md),
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md),
and
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
each do one job, and none of them need the others – call whichever one
you actually need. The first time you’re working with a new dataset,
though, run all three in order – diagnose, clean, extract – so you
actually see what’s wrong with the data before it quietly turns into
missing rows.

- **[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md)
  – diagnose.** Read-only. Scans the input and returns a `karyo_check`
  tibble flagging issue types (dirty markers, chimeric separators,
  structural errors), plus `fixable`/`unfixable` summary columns.
- **[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
  – clean.** The only function that rewrites a string. Fixes what’s
  fixable and returns a narrow `karyo_preprocessed` tibble (`original`,
  `preprocessed`, `status`).
- **[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
  – extract.** Feature extraction only – no cleaning of its own.
  `on_issues = "preprocess"` can invoke the same check/fix logic
  internally, so a single call is often enough; chain the three
  explicitly when you also want the diagnostic or cleaned tibble in its
  own right. See
  [`vignette("data-cleaning")`](https://lucalanino.github.io/karyoparser/articles/data-cleaning.md)
  for the full `on_issues` reference.

## Chaining pipeline steps

[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md)
and
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
tag their output with an extra S3 class, so you can pipe straight into
the next step without re-specifying arguments:

| Function | Class | Notes |
|----|----|----|
| [`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md) | `karyo_check` (+ `tbl_df`) | Accepted as [`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md) input; the assessment is reused instead of re-run. |
| [`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md) | `karyo_preprocessed` (+ `tbl_df`) | Accepted as [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md) input; the `preprocessed` column and cached id/issue info are reused automatically. |
| [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md) | plain `tbl_df` | No custom class. |

[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)’s
output also carries a `karyoparser_version` attribute
(`attr(result, "karyoparser_version")`) set to the installed package
version, so you can trace a saved result back to whatever version
produced it.

## A first pass, with `example_karyotypes`

`example_karyotypes` ships with the package – 102 synthetic ISCN strings
covering one example per `myeloid_rules` flag, general structural
aberrations, aneuploidy-only cases, balanced/unbalanced translocations,
complex and monosomal karyotypes, chimeric/multi-clone cases, and a
batch of messy or out-of-scope strings.

``` r

dplyr::glimpse(example_karyotypes)
#> Rows: 102
#> Columns: 2
#> $ sample_id <chr> "EX001", "EX002", "EX003", "EX004", "EX005", "EX006", "EX007…
#> $ karyotype <chr> "46,XX[20]", "46,XY[20]", "46,XY", "46,XX,t(15;17)(q24;q21)[…

ck <- check_karyo(
  example_karyotypes,
  karyotype_column = "karyotype",
  id_column = "sample_id"
)
dplyr::count(ck, fixable, unfixable)
#> # A tibble: 3 × 3
#>   fixable unfixable     n
#>     <int>     <int> <int>
#> 1       0         0    89
#> 2       0         1     4
#> 3       1         0     9
```

Fixable rows are formatting artifacts (stray whitespace, HTML entities,
a missing comma, …); unfixable rows are structural problems (no
chromosome count, mosaicism, …) that always come out `NA`.
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
only touches the fixable kind. See
[`vignette("data-cleaning")`](https://lucalanino.github.io/karyoparser/articles/data-cleaning.md)
for the full breakdown of issue types and how `on_issues` handles each.

``` r

pp <- preprocess_karyo(ck)
result <- parse_karyo(pp, on_issues = "preprocess")
#> Chimeric: 3 (on_chimeric = "default").
result |>
  dplyr::select(sample_id, chromosome_count, general_translocation, complex_karyotype)
#> # A tibble: 102 × 4
#>    sample_id chromosome_count general_translocation complex_karyotype
#>    <chr>                <int>                 <int>             <int>
#>  1 EX001                   46                     0                 0
#>  2 EX002                   46                     0                 0
#>  3 EX003                   46                     0                 0
#>  4 EX004                   46                     1                 0
#>  5 EX005                   46                     1                 0
#>  6 EX006                   46                     0                 0
#>  7 EX007                   46                     1                 0
#>  8 EX008                   46                     1                 0
#>  9 EX009                   46                     1                 0
#> 10 EX010                   46                     0                 0
#> # ℹ 92 more rows
```

Naming the columns once at
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md)
is enough –
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
and
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
inherit them from the object they are handed. Passing the data frame
straight to `parse_karyo(..., on_issues = "preprocess")` does the same
check-and-fix work in one call, with the columns named there instead.

## Data frame input

[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
also takes a data frame directly. Name the karyotype column with
`karyotype_column`; columns are never guessed, so that argument is
required. `id_column` is optional – give it to carry an identifier into
the output, leave it out if the data has none.

``` r

df <- tibble::tibble(
  sample_id = c("S1", "S2", "S3"),
  iscn = c("46,XX", "47,XY,+21[10]/46,XY[5]", "46,XX,t(9;22)(q34;q11)[20]")
)

parse_karyo(df)
#> Error:
#> ! `karyotype_column` must be given when parse_karyo() is passed a data frame. 
#> Available columns: "sample_id", "iscn". 
#> e.g. parse_karyo(data, karyotype_column = "<col>")
```

The error lists the available columns, so the fix is usually a
copy-paste:

``` r

result <- parse_karyo(
  df,
  karyotype_column = "iscn",
  id_column = "sample_id",
  verbose = FALSE
)
result |>
  dplyr::select(sample_id, original_karyotype, chromosome_count, tris21)
#> # A tibble: 3 × 4
#>   sample_id original_karyotype         chromosome_count tris21
#>   <chr>     <chr>                                 <int>  <int>
#> 1 S1        46,XX                                    46      0
#> 2 S2        47,XY,+21[10]/46,XY[5]                   47      1
#> 3 S3        46,XX,t(9;22)(q34;q11)[20]               46      0
```

If the id column has duplicated values – serial karyotypes for one
patient, say – you get a warning rather than an error. Parsing is
unaffected, because rows are matched on the karyotype string and never
on the id, but joins on that column downstream may fan out.

## Output reference

### Output columns

Each row is one input karyotype; columns fall into a handful of groups
by provenance. The full output is wide – 170 columns for the default
`myeloid_rules` – so pass `columns` to
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
to keep just the groups you want (see
[`vignette("custom-rules")`](https://lucalanino.github.io/karyoparser/articles/custom-rules.md)
for details). Metadata and status columns are always included
regardless.

| Group | Columns | Type | Description |
|----|----|----|----|
| Metadata | `original_karyotype`, `preprocessed_karyotype` | character | Raw input; cleaned/normalized string that was parsed (NA for unfixable rows) |
| Classification | `normal_karyotype` | integer 0/1 | 1 if `46,XX` or `46,XY` exactly |
| Classification | `complex_karyotype` | integer 0/1 | 1 if \>= 3 distinct aberrations across clones |
| Classification | `monosomal_karyotype` | integer 0/1 | 1 if \>= 2 autosomal monosomies, or \>= 1 autosomal monosomy + \>= 1 structural aberration (any `general_*` flag except `general_marker` and `general_dmin`). Sex-chromosome monosomies never count. Forced to 0 for CBF-AML rows – `t(8;21)(q22;q22)`, `inv(16)(p13q22)`, `t(16;16)(p13;q22)` – even when the criteria are met |
| Rule flags | *(one per `myeloid_rules` entry)* | integer 0/1 | Specific lesions – see Aberration flags below |
| General flags | `general_translocation`, `general_deletion`, `general_inversion`, `general_addition`, `general_dicentric`, `general_isodicentric`, `general_pseudodicentric`, `general_isochromosome`, `general_ring`, `general_insertion`, `general_duplication`, `general_triplication`, `general_marker`, `general_dmin`, `general_derivative`, `balanced_translocation`, `unbalanced_translocation` | integer 0/1 | Universal structural-aberration detections, incl. translocation balance |
| Aneuploidy | `mono1`–`mono22`, `monoX`, `monoY`; `tris1`–`tris22`, `trisX`, `trisY` | integer 0/1 | Whole-chromosome loss/gain from `-`/`+` tokens |
| Summary | `comma_count_aberrations` | integer | Aberration count (max across clones; idem-expanded) |
| Summary | `distinct_aberrations` | integer | Distinct aberration tokens pooled across clones; the count `complex_karyotype` thresholds at 3 |
| Summary | `chromosome_count` | integer | Count from the most abnormal clone that clears `min_metaphases` (see below); NA if none does |
| Summary | `total_metaphases` | integer | Sum of bracket counts; NA if no brackets |
| Derived loss | `unbal_partial_loss_<arm>` (one per chromosome arm), `unbal_partial_loss` | integer 0/1 | Partial arm loss implied by an unbalanced der translocation |
| Status | `fixable_error`, `unfixable_error` | integer 0/1 | Row had a fixable / unfixable issue (unfixable rows are all NA) |
| Status | `chimeric_karyotype` | integer 0/1 | 1 if input contained a `//` chimeric separator |
| Status | `chimeric_clone` | character | Which clone was parsed for a chimeric row (`"host"`, `"donor"`, or NA) – see [`vignette("chimeric-karyotypes")`](https://lucalanino.github.io/karyoparser/articles/chimeric-karyotypes.md) |

### Which clone `chromosome_count` comes from

A karyotype can describe several clones, so
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
has to pick one to report a count for. It takes the **most abnormal**
clone – the one whose count is furthest from 46 – among those large
enough to be worth trusting.

“Large enough” is the `min_metaphases` argument, and you can change it:

``` r

k <- "25,X[3]/46,XX[10]"

# Default: a clone needs >= 2 metaphases, so the 3-metaphase clone counts.
parse_karyo(k, verbose = FALSE)$chromosome_count
#> [1] 25

# Raise the bar and the small clone is ignored.
parse_karyo(k, verbose = FALSE, min_metaphases = 5)$chromosome_count
#> [1] 46
```

Two things to know:

- If **no** clone clears the threshold, `chromosome_count` is `NA` –
  even though `total_metaphases` still reports the brackets it found.
  Seeing `NA` next to a non-`NA` `total_metaphases` means “every clone
  was too small”, not “no count present”.
- The threshold only applies to rows where at least one clone carries a
  metaphase bracket. A string with no brackets anywhere
  (e.g. `46,XX/45,XY,-7`) uses all of its clones regardless of
  `min_metaphases`.

Set `min_metaphases = 0` to consider every clone that has a count at
all.

### Aberration flags

Every flag is `0` or `1`, and **every flag fires independently** –
there’s no priority or competition between them, so one token can light
up several flags at once. They come from three sources:

- **Specific rule flags** – one per `myeloid_rules` entry, fired when
  its regex matches a token (e.g. `t(9;22)(q34;q11)`, `del(5q)`).
- **General structural flags** (`general_*`) – universal,
  disease-agnostic detections (any translocation, deletion, inversion,
  addition, dicentric, isodicentric, pseudodicentric, isochromosome,
  ring, insertion, duplication, triplication, marker, double minutes,
  derivative). Always on, regardless of which rule set you pass.
- **Aneuploidy flags** (`mono*`/`tris*`) – whole-chromosome loss/gain
  from `-`/`+` tokens.

Since flags don’t compete, a token can trip **all** its matches at once.
`t(9;11)(p21;q23)` sets the specific `t_9_11_p21_q23`, the family flag
`t_v_11q23`, and `general_translocation`; `del(5q)` sets both `del_5q`
and `general_deletion`.

- **Family detectors** (`t_v_11q23`, `t_v_11p15`, `t_3q26_v`, `t_5q`,
  `t_12p`) deliberately co-fire with the more specific flags – they mean
  “any rearrangement involving this region/arm”.
- **`*_other` variants** (e.g. `t_9_22_other`) are the complement of the
  canonical breakpoints: they fire for the same chromosome pair at
  *non-canonical* breakpoints only – the exclusion lives in the regex
  itself, so a canonical token sets just the specific flag.

``` r

myeloid_rules
#> # A tibble: 43 × 2
#>    flag_name         regex                                                      
#>  * <chr>             <chr>                                                      
#>  1 t(15;17)(q24;q21) "t\\(15;17\\)\\((q24|q22);q21\\)|t\\(17;15\\)\\(q21;(q24|q…
#>  2 t(8;21)(q22;q22)  "t\\(8;21\\)\\((q22|q21(\\.3)?);q22\\)|t\\(21;8\\)\\(q22;(…
#>  3 inv(16)(p13q22)   "inv\\(16\\)\\(p13q22\\)"                                  
#>  4 t(16;16)(p13;q22) "t\\(16;16\\)\\(p13;q22\\)"                                
#>  5 t(9;11)(p21;q23)  "t\\(9;11\\)\\(p21;q23\\)|t\\(11;9\\)\\(q23;p21\\)"        
#>  6 t(6;9)(p22;q34)   "t\\(6;9\\)\\(p22;q34\\)|t\\(9;6\\)\\(q34;p22\\)"          
#>  7 inv(3)(q21q26)    "inv\\(3\\)\\(q21q26\\)"                                   
#>  8 t(3;3)(q21;q26)   "t\\(3;3\\)\\(q21;q26\\)"                                  
#>  9 t(9;22)(q34;q11)  "t\\(9;22\\)\\(q34;q11\\)|t\\(22;9\\)\\(q11;q34\\)"        
#> 10 t(1;3)(p36;q21)   "t\\(1;3\\)\\(p36;q21\\)|t\\(3;1\\)\\(q21;p36\\)"          
#> # ℹ 33 more rows
```

## Next steps

- [`vignette("data-cleaning")`](https://lucalanino.github.io/karyoparser/articles/data-cleaning.md)
  – the `on_issues` modes, the fixable/unfixable issue catalog, and
  [`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)’s
  cleaning pipeline in detail.
- [`vignette("chimeric-karyotypes")`](https://lucalanino.github.io/karyoparser/articles/chimeric-karyotypes.md)
  – how `//`-separated chimeric karyotypes are handled.
- [`vignette("risk-stratification")`](https://lucalanino.github.io/karyoparser/articles/risk-stratification.md)
  – IPSS-R and ELN 2022 cytogenetic risk categories via
  [`assign_risk()`](https://lucalanino.github.io/karyoparser/reference/assign_risk.md).
- [`vignette("custom-rules")`](https://lucalanino.github.io/karyoparser/articles/custom-rules.md)
  – writing your own rule tables and restricting
  [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)’s
  output with `columns`.
