# TODO

## Parsing & data quality

- **[P3] Match-time tolerance not applied everywhere**: `rule_match_text()` drives rules and `general_*` only. `distinct_aberrations` still counts `inv(3)(q21;q26)` and `inv(3)(q21q26)` as different, and translocation arm derivation ignores `?q21` breakpoints.
