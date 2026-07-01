# TODO

## Before 1.0

- **[P1] R/ reorganization**: One file per export (roxygen on top); split helpers.R into cohesive groups (`input.R`, `tokens.R`, etc.); consolidate shared check/preprocess internals; design with FISH reuse in mind (expand, don't rewrite)
- **[P2] `lymphoid_rules`**: Add exported `karyo_rules` constant for lymphoid neoplasms (a curated table of lymphoid-specific lesions, mirroring `myeloid_rules`). General structural detections are already universal, so no shared base table is needed.
- **[P2] Myeloid rule completeness**: Final pass over `myeloid_rules` 
- **[P3] CI / GitHub Actions**: `R-CMD-check` and `format-check` workflows exist. Remaining: decide post-1.0 branch strategy (`main` public after 1.0 -- which branches trigger checks?) and add a coverage workflow (`covr::codecov()`).
- **[P3] Codecov badge**: when repo goes public, add `CODECOV_TOKEN` secret and wire up codecov upload in `test-coverage.yaml`; add badge to README.

## After 1.0

- **pkgdown site**: `usethis::use_pkgdown()` + GitHub Pages. Write vignettes. Audit `.Rbuildignore`/`.gitignore` before building -- `pkgdown::build_site()` can leak ignored files into `docs/`.
- **FISH**: develop functions to check, preprocess and parse FISH strings.
