
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

result[, c("original_karyotype", "chromosome_count", "tris21", "t_9_22_q34_q11", "normal_karyotype")]
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

### The three pipeline functions

`check_karyo()`, `preprocess_karyo()`, and `parse_karyo()` each do one
job, and none of them need the others – call whichever one you actually
need. The first time you’re working with a new dataset, though, run all
three in order – diagnose, clean, extract – so you actually see what’s
wrong with the data before it quietly turns into missing rows.

- **`check_karyo()` – diagnose.** Read-only. Scans the input and returns
  a `karyo_check` tibble flagging issue types (dirty markers, chimeric
  separators, structural errors), plus `fixable`/`unfixable` summary
  columns.
- **`preprocess_karyo()` – clean.** The only function that rewrites a
  string. Fixes what’s fixable and returns a narrow `karyo_preprocessed`
  tibble (`original`, `preprocessed`, `status`).
- **`parse_karyo()` – extract.** Feature extraction only – no cleaning
  of its own. `on_issues = "preprocess"` can invoke the same check/fix
  logic internally, so a single call is often enough; chain the three
  explicitly when you also want the diagnostic or cleaned tibble in its
  own right.

Each function takes the previous step’s output directly – `karyo_check`
into `preprocess_karyo()`, `karyo_preprocessed` into `parse_karyo()` –
and reuses the cached assessment and id column instead of re-scanning.
More on that in [Chaining pipeline steps](#chaining-pipeline-steps).

### Chaining pipeline steps

`check_karyo()` and `preprocess_karyo()` tag their output with an extra
S3 class, so you can pipe straight into the next step without
re-specifying arguments:

| Function | Class | Notes |
|----|----|----|
| `check_karyo()` | `karyo_check` (+ `tbl_df`) | Accepted as `preprocess_karyo()` input; the assessment is reused instead of re-run. |
| `preprocess_karyo()` | `karyo_preprocessed` (+ `tbl_df`) | Accepted as `parse_karyo()` input; the `preprocessed` column and cached id/issue info are reused automatically. |
| `parse_karyo()` | plain `tbl_df` | No custom class. |

`parse_karyo()`’s output also carries a `karyoparser_version` attribute
(`attr(result, "karyoparser_version")`) set to the installed package
version, so you can trace a saved result back to whatever version
produced it.

### Data frame input

`parse_karyo()` also takes a data frame directly. It auto-detects the
karyotype column (common names like `karyotype`, `iscn`, `karyo`, …) and
an optional id column (`sample_id`, `patient_id`, `id`, `mrn`, …).

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

### Fixable vs. unfixable

Every issue `check_karyo()` finds is either **fixable** (a formatting
artifact – stray whitespace, HTML entities, a missing comma, a FISH
suffix, …) or **unfixable** (a structural problem, like no chromosome
count). `preprocess_karyo()` only touches the fixable kind. Unfixable
rows always come out `NA`, no matter which `on_issues` mode you use –
that’s usually where an unexpected wall of `NA`s in your output comes
from. Check `check_karyo()`’s `fixable`/`unfixable` columns before you
parse, so you know which rows will vanish and why.

## Usage

### Example Dataset

`example_karyotypes` ships with the package – 100 synthetic ISCN strings
so you can try everything above without bringing your own data. It
covers one example per `myeloid_rules` flag, general structural
aberrations, aneuploidy-only cases, balanced/unbalanced translocations,
complex and monosomal karyotypes, plus a batch of messy or out-of-scope
strings (leading dots, HTML entities, FISH suffixes,
`mos`/`ncSCA`/constitutional markers, …) for exercising `check_karyo()`
and `preprocess_karyo()`.

``` r
example_karyotypes
#> # A tibble: 100 × 2
#>    sample_id karyotype                  
#>    <chr>     <chr>                      
#>  1 EX001     46,XX[20]                  
#>  2 EX002     46,XY[20]                  
#>  3 EX003     46,XY                      
#>  4 EX004     46,XX,t(15;17)(q24;q21)[20]
#>  5 EX005     46,XY,t(8;21)(q22;q22)[18] 
#>  6 EX006     46,XX,inv(16)(p13q22)[15]  
#>  7 EX007     46,XY,t(16;16)(p13;q22)[12]
#>  8 EX008     46,XX,t(9;11)(p21;q23)[14] 
#>  9 EX009     46,XY,t(6;9)(p22;q34)[16]  
#> 10 EX010     46,XX,inv(3)(q21q26)[10]   
#> # ℹ 90 more rows

ck <- check_karyo(example_karyotypes)
table(fixable = ck$fixable, unfixable = ck$unfixable)
#>        unfixable
#> fixable  0  1
#>       0 87  5
#>       1  8  0

result <- parse_karyo(example_karyotypes, on_issues = "preprocess", verbose = FALSE)
#> Chimeric: 3 (on_chimeric = "default").
result[, c("sample_id", "chromosome_count", "general_translocation", "complex_karyotype")]
#> # A tibble: 100 × 4
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
#> # ℹ 90 more rows
```

## Output Reference

### Output Columns

Each row is one input karyotype; columns fall into a handful of groups
by provenance. The full output is wide – 167 columns for the default
`myeloid_rules` – so pass `columns` to `parse_karyo()` to keep just the
groups you want, e.g. `columns = c("classification", "general")`.
Metadata and status columns are always included regardless. Valid
values: `"classification"`, `"rule"`, `"general"`, `"aneuploidy"`,
`"summary"`, `"der_loss"` (matching the groups below).

| Group | Columns | Type | Description |
|----|----|----|----|
| Metadata | `original_karyotype`, `preprocessed_karyotype` | character | Raw input; cleaned/normalized string that was parsed (NA for unfixable rows) |
| Classification | `normal_karyotype` | integer 0/1 | 1 if `46,XX` or `46,XY` exactly |
| Classification | `complex_karyotype` | integer 0/1 | 1 if \>= 3 distinct aberrations across clones |
| Classification | `monosomal_karyotype` | integer 0/1 | 1 if \>= 2 autosomal monosomies, or \>= 1 monosomy + \>= 1 structural aberration |
| Rule flags | *(one per `myeloid_rules` entry)* | integer 0/1 | Specific lesions – see [Aberration Flags](#aberration-flags) |
| General flags | `general_translocation`, `general_deletion`, `general_inversion`, `general_addition`, `general_dicentric`, `general_isodicentric`, `general_pseudodicentric`, `general_isochromosome`, `general_ring`, `general_insertion`, `general_duplication`, `general_triplication`, `general_marker`, `general_derivative`, `balanced_translocation`, `unbalanced_translocation` | integer 0/1 | Universal structural-aberration detections, incl. translocation balance |
| Aneuploidy | `mono1`–`mono22`, `monoX`, `monoY`; `tris1`–`tris22`, `trisX`, `trisY` | integer 0/1 | Whole-chromosome loss/gain from `-`/`+` tokens |
| Summary | `comma_count_aberrations` | integer | Aberration count (max across clones; idem-expanded) |
| Summary | `chromosome_count` | integer | Count from the most abnormal eligible clone |
| Summary | `total_metaphases` | integer | Sum of bracket counts; NA if no brackets |
| Derived loss | `unbal_partial_loss_<arm>` (one per chromosome arm), `unbal_partial_loss` | integer 0/1 | Partial arm loss implied by an unbalanced der translocation |
| Status | `fixable_error`, `unfixable_error` | integer 0/1 | Row had a fixable / unfixable issue (unfixable rows are all NA) |
| Status | `chimeric_karyotype` | integer 0/1 | 1 if input contained a `//` chimeric separator |
| Status | `chimeric_clone` | character | Which clone was parsed for a chimeric row (`"host"`, `"donor"`, or NA) |

### Aberration Flags

Every flag is `0` or `1`, and **every flag fires independently** –
there’s no priority or competition between them, so one token can light
up several flags at once. They come from three sources:

- **Specific rule flags** – one per `myeloid_rules` entry, fired when
  its regex matches a token (e.g. `t(9;22)(q34;q11)`, `del(5q)`).
- **General structural flags** (`general_*`) – universal,
  disease-agnostic detections (any translocation, deletion, inversion,
  addition, dicentric, isodicentric, pseudodicentric, isochromosome,
  ring, insertion, duplication, triplication, marker, derivative).
  Always on, regardless of which rule set you pass.
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

Here’s the full default rule set:

``` r
myeloid_rules
#> # A tibble: 42 × 2
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
#> # ℹ 32 more rows
```

## Custom Rules

``` r
my_rules <- validate_rules(dplyr::bind_rows(
  as.data.frame(myeloid_rules),
  data.frame(
    flag_name = "t(X;18)(p11;q11)",
    regex = "t\\(X;18\\)\\(p11;q11\\)|t\\(18;X\\)\\(q11;p11\\)"
  )
))

result <- parse_karyo(df, rules = my_rules)
```

Rule columns: `flag_name` and `regex`. Rules fire independently when
their `regex` matches a token, with no priority ranking – so if you need
a flag to *not* fire in some case (e.g. a “non-canonical breakpoints
only” variant), encode that exclusion in the `regex` itself (a negative
lookahead works well). `flag_name` must be unique – `validate_rules()`
errors on duplicates, whether within one table or across tables combined
via `list()` – so alternative patterns for the same flag belong in one
`regex` (joined with `|`), not separate rows. General structural
categories (translocations, deletions, dicentrics, …) are already
detected universally by the parser, so you don’t need rules for those.

`flag_name` can use ISCN-style punctuation for readability – the output
column is a sanitized version, with runs of non-alphanumeric characters
replaced by `_` (`t(X;18)(p11;q11)` becomes `t_X_18_p11_q11`). That
keeps output columns free of characters (`(`, `)`, `;`) that would
otherwise force backtick-quoting (`` r$`t(X;18)(p11;q11)` ``) in
interactive use. `validate_rules()` also catches the rare case where two
different `flag_name` values sanitize down to the same column name.

## What Is Not Handled

- **Sub-band breakpoints**: bands are stripped before matching, so
  precise breakpoints aren’t stored.
- **Copy number \> 1**: gain/loss is binary – `+8,+8` still just gives
  `tris8 = 1`.
- **Partial gain, and offsetting partial loss**: an unbalanced
  `der(a)t(a;b)` records the implied partial loss
  (`unbal_partial_loss_*`), but never the matching gain – there’s no
  `unbal_partial_gain_*`. The loss comes straight from the der’s own
  breakpoints, not from reconciling the whole karyotype, so a
  co-occurring gain that would offset it (e.g. `+9` alongside
  `der(9)t(9;22)(q34;q11)`) doesn’t suppress the loss flag –
  `unbal_partial_loss_9q` still fires even though distal 9q is really at
  two copies.
- **Sex chromosome syndromes**: we track plain monosomy/trisomy, not
  syndromes – no dedicated Turner or Klinefelter flag.
- **`idem` across clones**: `idem` gets expanded for aberration
  counting, but flag and balance checks stay clone-local. So a
  reciprocal der pair split across an `idem` boundary (one der in the
  stemline, its partner in an `idem` subclone) reads as
  `unbalanced_translocation` only – we miss the `balanced_translocation`
  the pair actually implies.
- **Three-way translocations**: a complete `t(a;b;c)` only trips the
  generic `general_translocation` flag and is assumed balanced – no
  specific-breakpoint rule fires for it. A lone `der()` from a three-way
  still gets classified balanced/unbalanced, but we don’t derive partial
  loss from it; that’s limited to plain two-partner `der(a)t(a;b)`.
- **Non-myeloid panels**: the bundled rules cover AML/MDS/MPN/CML.
  Lymphoid and solid-tumor lesions aren’t in scope.
- **Mosaicism (`mos`)**: flagged `mosaic_karyotype` and left unparsed.
- **Non-clonal single-cell abnormalities (`ncSCA`)**: flagged
  `non_clonal_sca` and left unparsed.
- **Constitutional abnormalities (`c`)**: a constitutional suffix like
  `47,XXYc` gets flagged `constitutional_sex_complement` instead of
  being parsed as an acquired change.
- **Array CGH / SNP arrays**: `seq[GRCh38]`-style notation isn’t parsed.
- **FISH-only results**: standalone FISH results outside an ISCN
  karyotype string aren’t processed.
