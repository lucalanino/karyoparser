# Check Karyotype Strings for Issues

Checks raw karyotype strings for formatting artifacts and structural
errors. Read-only: this is a diagnostic step and never modifies a
karyotype string (cleaning happens in
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md)).
Returns a wide-format tibble with one row per input string: each
unfixable issue type is a column (`0`/`1`), plus summary `fixable` and
`unfixable` columns. Fixable issue types are not broken out into columns
since
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
under `on_issues = "preprocess"` resolves them automatically; the
`fixable` column flags rows affected by any of them. When
`verbose = TRUE`, prints a count summary and a per-unfixable-issue-type
breakdown.

## Usage

``` r
check_karyo(
  karyotypes,
  karyotype_column = NULL,
  id_column = NULL,
  verbose = FALSE
)
```

## Arguments

- karyotypes:

  Character vector of karyotype strings, or a data frame containing a
  karyotype column. For a data frame, name the karyotype column with
  `karyotype_column`; optionally name an id column with `id_column`,
  which is then included as the first column of the output and
  propagated through subsequent pipeline steps.

- karyotype_column:

  Character. Name of the karyotype column. **Required** when
  `karyotypes` is a data frame – columns are never guessed. Ignored when
  `karyotypes` is a character vector.

- id_column:

  Character. Name of the id column. Optional: `NULL` (default) means the
  data has no identifier, and one is never inferred. When given it is
  placed first in the output and propagated downstream. Duplicate ids
  are warned about, not rejected: rows are matched on the karyotype
  string and never on the id, so parsing is unaffected – but downstream
  joins on the column may fan out.

- verbose:

  Logical. If `TRUE`, prints a count summary and a
  per-unfixable-issue-type breakdown. Default `FALSE`.

## Value

A `karyo_check` tibble with `length(karyotypes)` rows (or
`nrow(karyotypes)` when input is a data frame). Columns: optional id
column (first, when given), `karyotype` (full input string), `fixable`,
`unfixable`, then one integer column per unfixable issue type
(alphabetical). The tibble can be passed directly to
[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md),
which will reuse the cached assessment and propagate the id column
without re-scanning.

## Details

Fixable issue types (resolvable by
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
under `on_issues = "preprocess"`): `unicode_notation`,
`non_ascii_homoglyph`, `fullwidth_punctuation`, `embedded_newline`,
`html_entities`, `leading_dot`, `count_sex_separator`, `fish_notation`,
`trailing_narrative`, `midstring_linewrap`, `missing_sex_comma`,
`mar_space`, `non_clonal_sca`, `chimeric_separator`,
`zero_host_chimera`.

Unfixable issue types (always returned as NA by
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)):
`multiple_chimeric_separator`, `empty`, `no_chromosome_count`,
`updated_iscn`, `unbalanced_parentheses`, `unbalanced_brackets`,
`no_sex_complement`, `single_token`, `constitutional_sex_complement`,
`mosaic_karyotype`, `invalid_idem`, `unparseable_bracket`,
`stray_non_ascii`.

See
[`vignette("data-cleaning")`](https://lucalanino.github.io/karyoparser/articles/data-cleaning.md)
for what each issue type means, worked examples, and how to resolve the
unfixable ones.

## See also

[`preprocess_karyo()`](https://lucalanino.github.io/karyoparser/reference/preprocess_karyo.md),
[`vignette("data-cleaning")`](https://lucalanino.github.io/karyoparser/articles/data-cleaning.md)
