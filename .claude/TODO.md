# TODO

## Before 1.0

- **[P1] R/ reorganization**: One file per export (roxygen on top); split helpers.R into cohesive groups (`input.R`, `tokens.R`, etc.); consolidate shared check/preprocess internals; design with FISH reuse in mind (expand, don't rewrite)
- **[P2] `lymphoid_rules`**: Add exported `karyo_rules` constant for lymphoid neoplasms 
- **[P2] Myeloid rule completeness**: Final pass over `myeloid_rules` 
- **[P2] Unbalanced partial-loss follow-ups**: Loss derivation is limited to simple single-junction `der(a)t(a;b)`; the remainders below are out-of-scope for now (rationale documented at `classify_translocation_balance()`/`derive_unbalanced_loss()` in `R/parse_karyo.R`):
  - (a) widen curated `unbal_loss_*` arm set (currently 5q/7q/11q/12p/13q/17p/20q) or expose all arms;
  - (b) copy-number-aware GAIN derivation (needs whole-karyotype reasoning);
  - (c) der of a three-way `t(a;b;c)` (e.g. `der(9)t(9;22;11)`) gets no balanced/unbalanced flag and no loss derivation -- needs a separate 3-partner code path. It still fires `derivative_chromosome` and counts toward `complex_karyotype`/`monosomal_karyotype`; a complete `t(a;b;c)` fires `general_translocation` (assumed balanced) but no specific-breakpoint flag;
  - (d) reciprocal ders split across an `idem` subclone report `unbalanced=1` only, missing the `balanced=1` the idem-inherited pair implies -- needs `idem` token-expansion.
- **[P3] CI / GitHub Actions**: `R-CMD-check` and `format-check` workflows exist. Remaining: decide post-1.0 branch strategy (`main` public after 1.0 -- which branches trigger checks?) and add a coverage workflow (`covr::codecov()`).
- **[P3] Codecov badge**: when repo goes public, add `CODECOV_TOKEN` secret and wire up codecov upload in `test-coverage.yaml`; add badge to README.
- **[P3] README**: `README.Rmd` is in place; flesh out remaining gaps -- notably guidance on building new rule tables.
- **[P3] Revisit zero_host_chimera classification under `on_chimeric="host"`**: currently treated as **unfixable** (status `"unfixable"`, `unfixable_error=1`, NA output) since the host-only policy yields no parseable clone. Alternative: keep it classified as a fixable issue (`fixable_error=1`) that simply produces NA output. Decide which is the more intuitive contract.

## After 1.0

- **pkgdown site**: `usethis::use_pkgdown()` + GitHub Pages. Write vignettes. Audit `.Rbuildignore`/`.gitignore` before building -- `pkgdown::build_site()` can leak ignored files into `docs/`.
- **FISH**: develop functions to check, preprocess and parse FISH strings.
