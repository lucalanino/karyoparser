# Default Parsing Rules for Myeloid Neoplasms

A `karyo_rules` tibble of regex-based rules for myeloid-specific
recurrent lesions: specific translocations/inversions, variable-partner
translocations, and chromosome-arm-specific deletions/additions. Passed
to
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
via the `rules` argument by default. Use
[`validate_rules()`](https://lucalanino.github.io/karyoparser/reference/validate_rules.md)
to build a custom rule set with the same schema.

## Usage

``` r
myeloid_rules
```

## Format

A `karyo_rules` tibble with columns: `flag_name`, `regex`.

## Details

General, disease-agnostic structural aberrations (dicentrics, rings,
insertions, generic translocations/deletions, derivatives, etc.) are NOT
rules: they are detected universally by the parser and surfaced as the
`general_*` output columns, independently of the rule set in use and of
any specific rule firing on the same token.

## See also

[`validate_rules()`](https://lucalanino.github.io/karyoparser/reference/validate_rules.md)
