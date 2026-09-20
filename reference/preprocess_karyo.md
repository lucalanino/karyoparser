# Clean and Normalize ISCN Karyotype Strings

The complete cleaning and normalization pipeline: fixes every fixable
issue type (see
[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md))
in a fixed order, then runs `normalize_iscn()` as a final pass. Output
is fully normalized and parser-ready. This is the only place in the
workflow where cleaning or normalization occurs –
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
does none. See
[`vignette("data-cleaning")`](https://lucalanino.github.io/karyoparser/articles/data-cleaning.md)
for the exact fix order and worked examples.

## Usage

``` r
preprocess_karyo(
  karyotypes,
  karyotype_column = NULL,
  id_column = NULL,
  on_chimeric = c("default", "host", "donor"),
  verbose = FALSE
)
```

## Arguments

- karyotypes:

  Character vector of raw karyotype strings, the `karyo_check` tibble
  returned by
  [`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md),
  or a data frame containing a karyotype column. When a `karyo_check`
  tibble is supplied, the assessment it already computed is reused
  directly – no re-scanning – and any id column named upstream is
  propagated automatically. When a plain data frame is supplied, name
  the karyotype column with `karyotype_column`.

- karyotype_column:

  Character. Name of the karyotype column. **Required** when
  `karyotypes` is a plain data frame – columns are never guessed.
  Ignored for character vector or `karyo_check` input.

- id_column:

  Character. Name of the id column. Optional: `NULL` (default) means the
  data has no identifier, and one is never inferred. Ignored for
  character vector or `karyo_check` input, where the id is propagated
  from the upstream object's attributes. Duplicate ids are warned about,
  not rejected: rows are matched on the karyotype string and never on
  the id, so parsing is unaffected – but downstream joins on the column
  may fan out.

- on_chimeric:

  Character string controlling which clone is kept for chimeric
  karyotypes (those containing a `//` separator). One of:

  - `"default"`: host clone (before `//`) for normal chimeras; donor
    clone (after `//`) for zero-host chimeras (`.//` prefix).

  - `"host"`: host clone only; zero-host chimeras become NA (the NA is
    policy-induced, so the row stays classified as fixable, not
    unfixable).

  - `"donor"`: everything after `//` (the donor population). Rows with
    two or more `//` separators are always unfixable. The choice is
    baked into the `preprocessed` column and recorded so
    [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
    can reuse it.

- verbose:

  Logical. If `TRUE`, prints a processing summary. Default `FALSE`.

## Value

A `karyo_preprocessed` tibble. Columns: optional id column (first, when
given or propagated), `original` (raw input, always unchanged),
`preprocessed` (fully normalized string, `NA_character_` for unfixable
rows), `status` (`"clean"`, `"fixed"`, or `"unfixable"`).

The returned tibble can be passed directly to
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
without specifying `karyotype_column` or `id_column` – both are carried
on the object's class and cached attributes.

## Details

`zero_host_chimera` strings (`.//` or `//` prefix) are donor-only
chimeras with no host metaphases. The leading `.//` prefix is stripped,
leaving the donor clone(s), which are then parsed (under
`on_chimeric = "host"` these rows are instead returned NA). A remaining
non-ASCII character with no known automatic fix (`stray_non_ascii`) is
detected but not modified; such rows are returned as `NA` (unfixable).

## See also

[`check_karyo()`](https://lucalanino.github.io/karyoparser/reference/check_karyo.md),
[`vignette("data-cleaning")`](https://lucalanino.github.io/karyoparser/articles/data-cleaning.md)
