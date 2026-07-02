# TODO

## Before 1.0

- **[P2] `lymphoid_rules`**: Add exported `karyo_rules` constant for lymphoid neoplasms (a curated table of lymphoid-specific lesions, mirroring `myeloid_rules`). General structural detections are already universal, so no shared base table is needed.
- **[P2] Myeloid rule completeness**: Final pass over `myeloid_rules` 
- **[P3] CI / GitHub Actions**: `R-CMD-check` and `format-check` workflows exist. Remaining: decide post-1.0 branch strategy (`main` public after 1.0 -- which branches trigger checks?).

## After 1.0

- **pkgdown site**: `usethis::use_pkgdown()` + GitHub Pages. Write vignettes. Audit `.Rbuildignore`/`.gitignore` before building -- `pkgdown::build_site()` can leak ignored files into `docs/`.
