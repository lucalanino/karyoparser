# TODO

## Parsing & data quality

- **[P1] `sl`/`sdl` not expanded in `comma_count_aberrations`**: `compute_comma_counts()` (`R/flags.R`) inherits the stemline count only for `idem`, so `44,sl,...` and `44,sdl1,...` clones undercount. Can misclassify `complex_karyotype` in borderline cases. `sdlN` must inherit from its parent sideline, not the stemline.
- **[P2] Ambiguous dirty patterns** (unfixable, need a decision): bracket sums (`[6+3]`: add the counts?), `idem*2` (doubled clone?), typos (`+marl1`, `+dert(`).
- **[P2] `//` before `idem`/`sl` treated as chimeric**: `...[7]//46,idem,...` is almost certainly a typo for `/`, but the `idem` clone is dropped as donor.
- **[P2] Gains/losses read per token, not relative to the parent clone**: in `45,XY,-1,...[11]/44,sl,+1,...` the `+1` restores the lost chromosome but sets `tris1`; same for `-add(12)` then `+12`. Design limit: needs clone-aware copy-number resolution. Scope before implementing.
- **[P2] `unbal_partial_loss` marks the partner's arm as lost**: `der(5)t(5;17)(q12;q12)` sets `unbal_partial_loss_17p`, but under ISCN a lone `der(a)t(a;b)` keeps both normal `b` homologs. Confirm intent in `derive_unbalanced_loss()`.
- **[P3] `psu dic(...)` also sets `general_dicentric`**: `general_dicentric` (`dic\\(`) matches inside `psu dic(`. Confirm whether co-firing is intended; if not, add a lookbehind.
- **[P3] Heteromorphisms counted as aberrations**: `9qh+` passes `unrecognized_token` but still counts toward `comma_count_aberrations`.
- **[P3] Match-time tolerance not applied everywhere**: `rule_match_text()` drives rules and `general_*` only. `distinct_aberrations` still counts `inv(3)(q21;q26)` and `inv(3)(q21q26)` as different, and translocation arm derivation ignores `?q21` breakpoints.
