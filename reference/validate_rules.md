# Validate and Create a Custom Rules Table

Validates a data frame against the required schema for karyotype parsing
rules and returns it as a `karyo_rules` object accepted by
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md).

## Usage

``` r
validate_rules(rules)
```

## Arguments

- rules:

  A data frame with columns: `flag_name` (character), `regex`
  (character). Each rule fires independently when its `regex` matches a
  token; there is no priority/competition between rules, so mutual
  exclusivity (where wanted) must be encoded in the `regex` itself.
  `flag_name` must be unique – each maps to exactly one output column,
  so alternative patterns for the same flag must be combined into a
  single `regex` (e.g. with `|`) rather than given as separate rows.
  `flag_name` can use ISCN-style punctuation (e.g. `"t(9;22)(q34;q11)"`)
  for readability – the output column name is a sanitized version, with
  runs of non-alphanumeric characters replaced by `_` (e.g.
  `t_9_22_q34_q11`). `flag_name` values that sanitize to the same column
  name are rejected.

## Value

A `karyo_rules` object (tibble subclass) accepted by
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md).
