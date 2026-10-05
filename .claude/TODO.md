# TODO

## Parsing & data quality

- **[P1] Capitalised indicators silently ignored**: `Del(5)(q15q33)` sets neither `del_5q` nor `general_deletion`, and no issue is reported. Add a fixable dirty pattern that lowercases known indicators (`.aberr_indicators_paren`) when followed by `(`.
- **[P1] `sl`/`sdl` not expanded in `comma_count_aberrations`**: `compute_comma_counts()` (`R/flags.R`) inherits the stemline count only for `idem`, so `44,sl,...` and `44,sdl1,...` clones undercount. Can misclassify `complex_karyotype` in borderline cases. `sdlN` must inherit from its parent sideline, not the stemline.
- **[P2] Gains/losses read per token, not relative to the parent clone**: in `45,XY,-1,...[11]/44,sl,+1,...` the `+1` restores the lost chromosome but sets `tris1`; same for `-add(12)` then `+12`. Design limit: needs clone-aware copy-number resolution. Scope before implementing.
- **[P3] `psu dic(...)` also sets `general_dicentric`**: `general_dicentric` (`dic\\(`) matches inside `psu dic(`. Confirm whether co-firing is intended; if not, add a lookbehind.
