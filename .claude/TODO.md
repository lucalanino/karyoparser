# TODO

## Before 1.0

- **[P1] tests**: create a helper test unit with the pk wrapper. Confirm then when running test(), the helper is loaded first
- **[P1] R/ reorganization**: One file per export (roxygen on top); split helpers.R into cohesive groups (`input.R`, `tokens.R`, etc.); consolidate shared check/preprocess internals; design with FISH reuse in mind (expand, don't rewrite)
- **[P2] `lymphoid_rules`**: Add exported `karyo_rules` constant for lymphoid neoplasms 
- **[P2] Myeloid rule completeness**: Final pass over `myeloid_rules` 
- **[P3] CI / GitHub Actions**: decide post-1.0 branch strategy first (`main` public after 1.0 — which branches trigger checks?). Bundle `R CMD check` + coverage (`covr::codecov()`) in one workflow.
- **[P3] Codecov badge**: when repo goes public, add `CODECOV_TOKEN` secret and wire up codecov upload in `test-coverage.yaml`; add badge to README.
- **[P3] README.md**: work on a public-ready readme, including info to build new rule tables

## After 1.0

- **pkgdown site**: `usethis::use_pkgdown()` + GitHub Pages. Convert README to `README.Rmd`, write vignettes. Audit `.Rbuildignore`/`.gitignore` before building — `pkgdown::build_site()` can leak ignored files into `docs/`.
- **FISH**: develop functions to check, preprocess and parse FISH strings.
