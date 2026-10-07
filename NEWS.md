# karyoparser 1.0.1

## Bug fixes

* Monosomy and trisomy flags (`mono*`/`tris*`) now fire only when the whole
  token is a gain or loss (`+8`, `-7`, `-Xx2`). Previously a leading copy count
  was read as a chromosome number, so `+1~3r` set `tris1` and `+2mar` or
  `+2dmin` set `tris2`. Constitutional (`+21c`) and uncertain (`+?8`) gains no
  longer set a trisomy flag.
* `general_ring` now detects unidentified rings with or without a count
  (`+r`, `+1~3r`), and `general_marker` detects numbered markers (`+mar1`).

## New features

* `check_karyo()` reports a new unfixable issue, `unrecognized_gain_loss`, for
  `+`/`-` tokens of no recognized shape (e.g. `+8q`, `+23`) in any clone, so
  they are surfaced instead of silently skipped. `parse_karyo()` returns `NA`
  for these rows.
* New fixable issue, `case_notation`: `preprocess_karyo()` now canonicalizes
  letter case, uppercasing sex chromosomes (`46,xy` -> `46,XY`, `-y`,
  `t(x;5)`) and lowercasing ISCN keywords (`Del(5)(Q13Q33)`, `MAR`, `[CP20]`).
  Lowercase karyotypes, previously unfixable and returned as `NA`, are now
  parsed; under the default `on_issues = "stop"`, `parse_karyo()` directs them
  to `preprocess_karyo()`. The lowercase `x` copy multiplier (`der(1)x2`) is
  left untouched.
