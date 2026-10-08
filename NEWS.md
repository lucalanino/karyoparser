# karyoparser 1.0.1

## Bug fixes

* Monosomy and trisomy flags (`mono*`/`tris*`) now fire only when the whole
  token is a gain or loss (`+8`, `-7`, `-Xx2`). Previously a leading copy count
  was read as a chromosome number, so `+1~3r` set `tris1` and `+2mar` or
  `+2dmin` set `tris2`. Constitutional (`+21c`) and uncertain (`+?8`) gains no
  longer set a trisomy flag.
* `general_ring` now detects unidentified rings with or without a count
  (`+r`, `+1~3r`), and `general_marker` detects numbered markers (`+mar1`).
* `inc` (incomplete karyotype) no longer counts toward
  `comma_count_aberrations` or `distinct_aberrations`.
* Rules and `general_*` patterns now fire on two common non-standard
  breakpoint forms: a `;` between the breakpoints of a single-chromosome
  rearrangement (`inv(16)(p13.1;q22)` sets `inv_16_p13q22`,
  `inv(3)(q21;q26)` sets `inv_3_q21q26`) and an uncertain breakpoint
  (`del(5)(?q13q31)` sets `del_5q`). Only matching is affected; the reported
  karyotype strings are unchanged.
* `unbal_partial_loss_*` no longer marks an arm of the translocation partner
  as lost. Per ISCN, a lone `der(a)t(a;b)` replaces one normal `a` while both
  normal `b` homologs remain, so only `a`'s arm beyond its breakpoint is lost:
  `der(5)t(5;17)(q12;q12)` now sets `unbal_partial_loss_5q` only, not
  `unbal_partial_loss_17p`. If you combined `unbal_partial_loss_17p` into
  abn(17p) for risk scoring, affected rows were wrongly scored adverse.
* `general_dicentric` no longer fires on `idic(...)` or `psu dic(...)`, which
  set only `general_isodicentric` and `general_pseudodicentric`. Each dicentric
  form now sets exactly one flag, as `general_isochromosome` already did.
* Spelling variants of one aberration (`del(5)(q13;q33)` and
  `del(5)(q13q33)`, or an uncertain `?q13`) now count once in
  `distinct_aberrations`, so they can no longer push a row into
  `complex_karyotype`. Translocation loss derivation also reads arms from
  uncertain breakpoints (`der(3)t(3;5)(?p13;q31)` sets
  `unbal_partial_loss_3p`).

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
* `check_karyo()` reports a new unfixable issue, `unrecognized_token`, for any
  aberration token with no recognized ISCN shape (free text, a bare number as
  in `46,XY,8,del(5)(q13q33)[20]`, an indicator with no breakpoints as in
  `47,XY,+8,del[20]`). Such tokens were previously
  counted as aberrations and could set `complex_karyotype`. Only token
  structure is checked, never band contents.
* `preprocess_karyo()` now strips a trailing angle-bracket annotation
  (`<AML>`) as `trailing_narrative`, and removes the dot of a `,.add(...)`
  token as `midstring_linewrap`.
* New fixable issue, `iso_indicator`: the non-ISCN `iso(` isochromosome
  indicator is rewritten to `i(` (`iso(17q)` -> `i(17q)`, setting `i_17q`).
* New fixable issues for wrong separators: `paren_comma` (a comma inside
  parentheses, `del(5)(q11,q33)`), `colon_separator` (`t(5:17)(q13:q11)`;
  groups in the detailed ISCN system are left alone) and `dot_separator`
  (`del(5)(q13q33).-7`, `45.idem`). These were previously split into
  fragments, read with too few flags, or silently dropped a token.
* More fixable issues: `underscore_prefix` (`46,_XY,_del(5)...`),
  `star_multiplier` (`+mar*2` -> `+marx2`) and `gain_loss_comma`
  (`-7+mar` -> `-7,+mar`). `paren_comma` reads a comma between band digits
  as a decimal comma (`q11,2` -> `q11.2`), and an unidentified marker, ring or
  double minute may now carry a copy multiplier (`+mar1x2`).
* `check_karyo()` gains an `unfixable_reason` column that names the
  offending token for each unfixable row (e.g. `unparseable_bracket: Cannot
  parse '[23+6]' as metaphase count`), so the reason is visible without
  inspecting internals.
* `check_karyo()` reports a new unfixable issue, `chimeric_no_count`, for a
  `//` that is not followed by a chromosome count (`46,XX[10]//XY[5]`, a
  trailing `//`). These rows were previously parsed as their host clone.
* `check_karyo()` reports a new unfixable issue, `invalid_breakpoint`, for a
  recognized token whose parentheses hold an impossible chromosome
  (`der(+)`) or breakpoint: an empty entry (`del(12)()`), a band without an
  arm (`t(5;17)(q12;21)`), or a breakpoint count that doesn't match the
  chromosomes (`t(5;2)(q11)`). ISCN uncertainty, `or` alternatives and
  uncertain ranges (`q13-14`) are accepted. `no_chromosome_count` now also
  fires when a later clone doesn't start with a chromosome count. These rows
  were previously parsed with partial or wrong flags.
* New fixable issue, `breakpoint_semicolon`: a missing `;` between two
  chromosomes' breakpoints is inserted (`t(9;16)(q34p13)` ->
  `t(9;16)(q34;p13)`).
