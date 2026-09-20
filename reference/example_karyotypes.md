# Synthetic Example Karyotypes

A synthetic dataset of 102 ISCN karyotype strings for trying out
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md),
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md),
and
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md).
Covers normal karyotypes, one example per `myeloid_rules` flag, general
structural aberrations with no matching specific rule, aneuploidy-only
cases, balanced/unbalanced translocations, complex and monosomal
karyotypes, multi-clone/chimeric/ploidy variety, and a range of messy or
out-of-scope strings (leading dots, HTML entities, FISH suffixes,
`mos`/`ncSCA`/ constitutional markers, ...) to demonstrate
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md)
and
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md).

## Usage

``` r
example_karyotypes
```

## Format

A tibble with 102 rows and columns:

- sample_id:

  character. Synthetic sample identifier (`"EX001"`, ...).

- karyotype:

  character. Raw ISCN karyotype string.

## See also

[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md),
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md),
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)
