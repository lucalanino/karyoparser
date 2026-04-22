# TODO

## Before 1.0

- **[P1] `lymphoid_rules`**: Add exported `karyo_rules` constant for lymphoid neoplasms 
- **[P1] Myeloid rule completeness**: Final pass over `myeloid_rules` 
- **[P2] CI / GitHub Actions**: decide post-1.0 branch strategy first (`main` public after 1.0 — which branches trigger checks?). Bundle `R CMD check` + coverage (`covr::codecov()`) in one workflow.
- **[P2] README.md**: work on a public-ready readme, including info to build new rule tables
## After 1.0

- **pkgdown site**: `usethis::use_pkgdown()` + GitHub Pages. Convert README to `README.Rmd`, write vignettes. Audit `.Rbuildignore`/`.gitignore` before building — `pkgdown::build_site()` can leak ignored files into `docs/`.
- **FISH**: develop functions to check, preprocess and parse FISH strings.
