# karyoparser <img src="man/figures/logo.png" align="right" height="138" alt="karyoparser logo" />

An R package for parsing ISCN karyotype strings into structured binary features. Designed for analysis of myeloid neoplasm-related chromosomal aberrations.

**Version**: 0.5.1

## Installation

```r
devtools::install_github("lucalanino/karyo-parser")
```

**Requires**: R >= 4.1.0

## Quick Start

```r
# 1. Install
devtools::install_github("lucalanino/karyo-parser")

# 2. Load
library(karyoparser)

# 3. Parse a character vector
result <- parse_karyo(
  c("46,XX", "47,XY,+21[10]/46,XY[5]", "46,XX,t(9;22)(q34;q11)[20]"),
  on_issues = "warn",
  verbose   = FALSE
)

# 4. Inspect key columns
result[, c("original_karyotype", "ploidy_category", "tris21", "t(9;22)(q34;q11)", "normal_karyotype")]
#> # A tibble: 3 × 5
#>   original_karyotype             ploidy_category tris21 `t(9;22)(q34;q11)` normal_karyotype
#>   <chr>                          <chr>            <int>              <int>            <int>
#> 1 46,XX                          diploid              0                  0                1
#> 2 47,XY,+21[10]/46,XY[5]         hyperdiploid         1                  0                0
#> 3 46,XX,t(9;22)(q34;q11)[20]     diploid              0                  1                0
```

That's it — one function call, one tibble. Each row is an input karyotype; each column is a binary flag (0/1), a ploidy category, or a count.

**Common next steps:**

```r
# Parse a data frame (karyotype column auto-detected)
df <- data.frame(
  sample_id = c("S1", "S2", "S3"),
  iscn      = c("46,XX", "47,XY,+21[10]/46,XY[5]", "46,XX,t(9;22)(q34;q11)[20]")
)
result <- parse_karyo(df)   # sample_id carried through automatically

# Dirty strings? Use on_issues = "fix" to auto-clean before parsing:
result <- parse_karyo(c(".46,XX", "47,XY,+21[10] .Lab note"), on_issues = "fix")

# Inspect what issues were found without parsing:
check_karyo(c(".46,XX", "bad string"))
#> # A tibble: 2 × 4
#>   row_index karyotype  issue_type   issue_detail
#>       <int> <chr>      <chr>        <chr>
#> 1         1 .46,XX     leading_dot  String starts with dot(s) before chromosome count
#> 2         2 bad string no_chrom...  ...

# Export
readr::write_csv(result, "parsed.csv")
```

## Usage

```r
library(karyoparser)

# Character vector
result <- parse_karyo(c("46,XX", "47,XY,+21", "46,XX,t(15;17)(q24;q21)"))

# Data frame — karyotype column auto-detected from common names (karyotype, iscn, karyo, ...)
result <- parse_karyo(df)

# Data frame with explicit column name
result <- parse_karyo(df, karyotype_column = "iscn")

# Silent, automated pipeline
result <- parse_karyo(df, verbose = FALSE)

# Export
readr::write_csv(result, "parsed.csv")
```

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `karyotypes` | — | Character vector or data frame |
| `karyotype_column` | `NULL` | Column name when input is a data frame (auto-detected if `NULL`) |
| `id_column` | `NULL` | ID column to carry through to output (auto-detected from sample_id, patient_id, id, mrn, etc.) |
| `rules` | `rules_table()` | Custom rules data frame (see [Custom Rules](#custom-rules)) |
| `.return` | `"tibble"` | Return type: `"tibble"` or `"data.frame"` |
| `verbose` | `TRUE` | Print validation reports and messages |
| `on_issues` | `"fix"` | `"fix"` (auto-clean dirty rows + truncate chimeric rows; print summary), `"warn"` (skip all issue rows → NA; print message), `"stop"` (error on any issue) |

## Output Columns

| Column | Type | Description |
|---|---|---|
| `original_karyotype` | character | Raw input string exactly as provided |
| `normalized_karyotype` | character | Input after preprocessing and ISCN normalization (what was actually parsed); NA for unfixable rows |
| `ploidy_category` | character | `diploid`, `hyperdiploid`, `high_hypodiploid`, `low_hypodiploid`, `near_haploid`, `other`, or `unknown` |
| `chromosome_count` | integer | Count from the most abnormal eligible clone |
| *(aberration flags)* | integer 0/1 | One column per rule — see [Aberration Flags](#aberration-flags) |
| `mono1`–`mono22`, `monoX`, `monoY` | integer 0/1 | Monosomy flags (from `−` tokens) |
| `tris1`–`tris22`, `trisX`, `trisY` | integer 0/1 | Trisomy flags (from `+` tokens) |
| `normal_karyotype` | integer 0/1 | 1 if `46,XX` or `46,XY` exactly |
| `total_metaphases` | integer | Sum of bracket counts; `NA` if no brackets |
| `comma_count_aberrations` | integer | Comma-separated aberrations (max across clones; idem-expanded) |
| `complex_karyotype` | integer 0/1 | 1 if ≥3 unique aberrations |
| `monosomal_karyotype` | integer 0/1 | 1 if ≥2 autosomal monosomies, or ≥1 monosomy + ≥1 structural aberration |
| `mixed_ploidy` | integer 0/1 | 1 if clones span different ploidy categories |
| `issues` | list | `NULL` for clean rows; tibble of `issue_type`/`issue_detail` for problem rows |

## How Rule Matching Works

Each karyotype string is split into comma-separated tokens (one per clone, idem-expanded). Rules are matched against each token independently.

Rules are organized into **competition groups** (translocation, inversion, deletion, addition, derivative, etc.). Within a group, only the highest-priority matching rule fires — so `t(9;22)(q34;q11)` (priority 100) suppresses `general_translocation` (priority 65) for the same token. Rules in **different groups co-fire independently**, so a complex token like `der(5)del(5)(q11q34)` sets both `del(5q) = 1` (deletion group) and `derivative_chromosome = 1` (derivative group). A token can therefore contribute to multiple output flags simultaneously.

Priority levels:

| Priority | Rule type |
|---|---|
| 100 | Specific translocations/inversions at canonical breakpoints |
| 95 | Variant-band translocations; isodicentric X |
| 90 | Variable-partner translocations |
| 85 | Chromosome-arm-specific aberrations; isodicentric; pseudodicentric |
| 80 | Dicentric |
| 65–70 | General structural (ring, insertion, duplication, triplication, general translocation) |
| 55–60 | Catch-all general rules (addition, inversion, deletion, marker, derivative) |

**CBF override**: Cases carrying `t(8;21)(q22;q22)`, `inv(16)(p13q22)`, or `t(16;16)(p13;q22)` are never classified as monosomal, per clinical guidelines, regardless of co-occurring monosomies.

## Aberration Flags

All flags output `0` or `1`. Monosomy/trisomy flags are derived from `−`/`+` tokens and listed separately under [Aneuploidy](#aneuploidy).

### Specific translocations (priority 100)

Both orientations are recognized (e.g. `t(15;17)` and `t(17;15)` both set the same flag). `counts_for_monosomal = FALSE` for CBF-AML lesions.

| Flag | Lesion |
|---|---|
| `t(15;17)(q24;q21)` | APL — *PML::RARA* |
| `t(8;21)(q22;q22)` | CBF-AML — *RUNX1::RUNX1T1* |
| `inv(16)(p13q22)` | CBF-AML — *CBFB::MYH11* |
| `t(16;16)(p13;q22)` | CBF-AML — *CBFB::MYH11* |
| `t(9;11)(p21;q23)` | *KMT2A::MLLT3* |
| `t(6;9)(p22;q34)` | *DEK::NUP214* |
| `inv(3)(q21q26)` | *GATA2::MECOM* |
| `t(3;3)(q21;q26)` | *GATA2::MECOM* |
| `t(9;22)(q34;q11)` | BCR-ABL1 |
| `t(1;3)(p36;q21)` | *RPN1::MECOM* |
| `t(1;22)(p13;q13)` | *RBM15::MRTFA* (acute megakaryoblastic) |
| `t(3;5)(q25;q35)` | *NPM1::MLF1* |
| `t(5;11)(q35;p15)` | *NUP98::NSD1* |
| `t(7;12)(q36;p13)` | *ETV6::MNX1* |
| `t(8;16)(p11;p13)` | *KAT6A::CREBBP* |
| `t(10;11)(p12;q14)` | *PICALM::MLLT10* |
| `t(11;12)(p15;p13)` | *NUP98::ETV6* |
| `t(16;21)(p11;q22)` | *FUS::ERG* |
| `t(16;21)(q24;q22)` | *RUNX1T3::RUNX1* |
| `inv(16)(p13q24)` | *CBFA2T3::GLIS2* (pediatric AML) |

### Variant-band translocations (priority 95)

Same chromosome pair as a canonical rule but at non-canonical sub-bands.

| Flag | Meaning |
|---|---|
| `t(6;9)_other` | t(6;9) at non-canonical bands |
| `t(9;11)_other` | t(9;11) at non-canonical bands |
| `t(9;22)_other` | t(9;22) at non-canonical bands |
| `inv(3)_other` | inv(3)(q?q?) not matching canonical q21q26 |
| `t(3;3)_other` | t(3;3) at non-canonical bands |
| `idic(X)(q13)` | Isodicentric X at q13 |

### Variable-partner translocations (priority 90)

| Flag | Meaning |
|---|---|
| `t(v;11p15)` | Any partner translocated to 11p15 (*NUP98* locus) |
| `t(v;11q23)` | Any partner translocated to 11q23 (*KMT2A* locus) |
| `t(3q26;v)` | Any partner translocated to 3q26 (*MECOM* locus) |

### Chromosome-arm-specific aberrations (priority 85)

| Flag | Meaning |
|---|---|
| `del(5q)` | Any deletion of chromosome 5 long arm |
| `t(5q)` | Any translocation involving 5q |
| `add(5q)` | Additional material of unknown origin on 5q |
| `del(7q)` | Any deletion of chromosome 7 long arm |
| `del(12p)` | Any deletion of chromosome 12 short arm |
| `t(12p)` | Any translocation involving 12p (*ETV6* locus) |
| `add(12p)` | Additional material of unknown origin on 12p |
| `del(13q)` | Any deletion of chromosome 13 long arm |
| `i(17q)` | Isochromosome 17q |
| `add(17p)` | Additional material of unknown origin on 17p |
| `del(17p)` | Any deletion of chromosome 17 short arm (*TP53* locus) |
| `del(20q)` | Any deletion of chromosome 20 long arm |
| `del(11q)` | Any deletion of chromosome 11 long arm |
| `isodicentric` | Any isodicentric chromosome (`idic(`) |
| `pseudodicentric` | Pseudodicentric chromosome (`psu dic(`) |

### General structural aberrations (priority 55–80)

Catch-all flags that fire when no more-specific rule matches within the same competition group.

| Flag | Matches | counts_for_monosomal |
|---|---|---|
| `dicentric` | `dic(` | yes |
| `ring_chromosome` | `r(` | yes |
| `insertion` | `ins(` | yes |
| `duplication` | `dup(` | yes |
| `triplication` | `trp(` | yes |
| `general_translocation` | `t(` not matched by a specific rule | yes |
| `general_inversion` | `inv(` not matched by a specific rule | yes |
| `general_deletion` | `del(` not matched by a specific rule | yes |
| `general_addition` | `add(` not matched by a specific rule | yes |
| `marker_chromosome` | `mar` | **no** |
| `derivative_chromosome` | `der(`, `ider(`, `+der(` | yes |

### Aneuploidy

Monosomy and trisomy flags (`mono1`–`mono22`, `monoX`, `monoY`, `tris1`–`tris22`, `trisX`, `trisY`) are detected from leading `−`/`+` tokens. Whole-chromosome gain/loss of all autosomes and sex chromosomes is captured. These are computed independently of the rule-matching pipeline and always present in the output.

## Preprocessing

Some data sources include dirty markers that the internal ISCN normalizer cannot fix: trailing clinical narrative, database prefixes, and HTML entities.

`parse_karyo()` detects these automatically via `check_karyo()`. The default `on_issues = "fix"` auto-cleans dirty rows and prints a summary of what was fixed. To clean manually before parsing, or to skip all issue rows instead:

```r
# Default: auto-fix dirty rows and truncate chimeric rows
result <- parse_karyo(raw_strings)

# Or clean manually, then parse:
clean  <- preprocess_karyo(raw_strings)
result <- parse_karyo(clean)

# Skip all issue rows (return NA) with a message:
result <- parse_karyo(raw_strings, on_issues = "warn")
```

Rules applied by `preprocess_karyo()` in order: trim whitespace → decode HTML entities (`&lt;`/`&gt;`/`&amp;`) → strip leading `.//` → strip leading dot before digit → strip trailing `] .text` narrative → collapse mid-string `, .token` line-wrap artifacts.

To inspect all issues (dirty markers + structural problems) without parsing:

```r
check_karyo(raw_strings)
# Returns a tibble with columns: row_index, karyotype, issue_type, issue_detail
```

Each output row includes an `issues` list-column: `NULL` for clean rows, or a tibble of issues for problem rows.

## Validation

The parser validates all input before parsing and reports issues by type:

| Issue type | Severity | Effect |
|---|---|---|
| `empty` | Error | Row gets NA values |
| `no_chromosome_count` | Error | Row gets NA values |
| `chimeric_separator` | Error | Row gets NA unless `on_issues = "fix"` (truncates to first clone) |
| `updated_iscn` | Error | Row always gets NA values |
| `unbalanced_parentheses` | Error | Row gets NA values |
| `unbalanced_brackets` | Error | Row gets NA values |
| `no_sex_complement` | Warning | Still parsed |
| `invalid_idem` | Warning | Still parsed (idem ignored) |
| `unparseable_bracket` | Warning | Bracket treated as 0 |

### `on_issues` behavior

| Value | Dirty rows | Chimeric rows (`//`) | Structural errors |
|---|---|---|---|
| `"fix"` (default) | `preprocess_karyo()` applied; still-dirty → NA | Truncated to first clone; still-invalid → NA | Always NA |
| `"warn"` | NA + message | NA + message | Always NA |
| `"stop"` | Error immediately | Error immediately | Error immediately |

Under `"fix"`, a message is printed summarising how many rows were cleaned and how many could not be fixed. Structural errors (unbalanced brackets, no chromosome count, `Updated ISCN` marker) are always NA regardless of `on_issues`; they cannot be auto-corrected.

## Custom Rules

```r
my_rules <- rules_table()

# Add a new rule
my_rules <- dplyr::bind_rows(my_rules, tibble::tibble(
  flag_name            = "t(X;18)(p11;q11)",
  regex                = "t\\(X;18\\)\\(p11;q11\\)|t\\(18;X\\)\\(q11;p11\\)",
  category             = "specific_tx",
  priority             = 100,
  counts_for_monosomal = TRUE,
  competition_group    = "translocation"
))

result <- parse_karyo(df, rules = my_rules)
```

Rule columns:

| Column | Description |
|---|---|
| `flag_name` | Output column name |
| `regex` | R-compatible regex matched against each token (bands stripped before matching) |
| `category` | Grouping label (`specific_tx`, `variable_partner`, `chromosome_specific`, `general`) |
| `priority` | Integer; higher = wins within its competition group when multiple rules match |
| `counts_for_monosomal` | Whether this aberration counts toward `monosomal_karyotype` |
| `competition_group` | Rules in the same group compete; rules in different groups co-fire independently |

## Special Notation

- **idem**: Refers to the stemline (clone 1). Subclones with `idem` inherit the stemline's aberration count for `comma_count_aberrations`. Excluded from unique aberration counts used for `complex_karyotype`.
- **sl**: Stemline marker. Excluded from unique aberration counts.
- **Composite karyotypes** (`cpN`): For ploidy classification, only clones with ≥5 metaphases are eligible. The most abnormal eligible clone (furthest from 46) determines the ploidy category.
- **Range notation** (`45~47,XX,...`): Treated as a single range karyotype; `normal_karyotype` is always 0.
- **Deduplication**: When the input contains duplicate karyotype strings, each unique string is parsed only once and results are broadcast back to all matching rows. No overhead when all strings are unique.

## What Is Not Handled

The following are outside the current scope. Some are deliberate design decisions; others are potential future extensions.

**Deliberate omissions:**

- **Sub-band breakpoint storage**: Cytogenetic bands are stripped before regex matching and not stored. The parser flags the *type* of aberration (e.g. `del(5q)`) but does not record precise breakpoints (e.g. `q13.1`).
- **Copy number beyond 0/1**: Gain and loss are binary. `+8,+8` yields `tris8 = 1`, not a copy-number count. Multiple copies of the same gain are not distinguished.
- **Sex chromosome numerical abnormalities** beyond simple monosomy/trisomy: 45,X is captured as `monoX = 1`; 47,XXY as `trisX = 1` and `trisY = 1`. No dedicated flags for Klinefelter, Turner, or sex chromosome polysomy syndromes.
- **Partial monosomy/trisomy via der/dup**: A `der(7)t(1;7)` carrying partial 7q loss is not auto-counted as a monosomy for chromosome 7.
- **Non-myeloid disease panels**: Rules are curated for myeloid neoplasms (AML, MDS, MPN, CML). Lymphoma-specific lesions (e.g. `t(14;18)`, `t(8;14)`) and solid-tumor translocations are absent.
- **Array CGH / SNP array notation**: The `seq[GRCh38]` format used in molecular karyotyping is not parsed.
- **FISH-only results**: Single-locus FISH findings expressed outside ISCN karyotype strings are not processed.

**Potential future extensions (not yet discussed):**

- *Per-clone output*: Splitting composite karyotypes into one row per clone, with a `clone_abundance` column (fraction of metaphases). Useful for clonal evolution analyses.
- *Cytogenetic risk group assignment*: Derived output columns for established schemas (ELN 2022, MRC, IPSS-R, IPSS-M) based on the parsed flags.
- *Gene-fusion annotation*: Mapping specific translocation flags to predicted gene partners (e.g. `t(9;22)` → BCR-ABL1), useful for linking cytogenetics to molecular data.
- *Lymphoid / solid-tumor rule sets*: Separate `rules_table()` variants for ALL, CLL, lymphoma, or sarcoma.
- *Gain/loss count columns*: Integer copy-number columns (e.g. `n_gains`, `n_losses`) distinct from the existing `comma_count_aberrations`.
- *Structural variant scoring*: A weighted aberration score integrating clinical weights (e.g. down-weighting sex chromosome loss in scoring models).
