# TODO

## Parsing & data quality

- **[P2] Gains/losses read per token, not relative to the parent clone**: in `45,XY,-1,...[11]/44,sl,+1,...` the `+1` restores the lost chromosome but sets `tris1`; same for `-add(12)` then `+12`. Design limit: needs clone-aware copy-number resolution. Scope before implementing.
- **[P2] `unbal_partial_loss` marks the partner's arm as lost**: `der(5)t(5;17)(q12;q12)` sets `unbal_partial_loss_17p`, but under ISCN a lone `der(a)t(a;b)` keeps both normal `b` homologs. Confirm intent in `derive_unbalanced_loss()`.
- **[P3] `psu dic(...)` also sets `general_dicentric`**: `general_dicentric` (`dic\\(`) matches inside `psu dic(`. Confirm whether co-firing is intended; if not, add a lookbehind.
- **[P3] Match-time tolerance not applied everywhere**: `rule_match_text()` drives rules and `general_*` only. `distinct_aberrations` still counts `inv(3)(q21;q26)` and `inv(3)(q21q26)` as different, and translocation arm derivation ignores `?q21` breakpoints.
