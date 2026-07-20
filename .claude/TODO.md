# TODO

## Before 1.0

- **[P2] Myeloid rule completeness**: Final pass over `myeloid_rules` 
- **[P3] Branch strategy (decided, not yet implemented)**: `main` stays the default branch (so untagged `pak::pak()`/`install_github()` installs always resolve to the last release -- confirmed via pak docs: "if `<detail>` is missing, the latest commit of the default branch is used"). Ongoing work moves to a `dev` branch (not set as default). Releases: PR `dev -> main`, require `R-CMD-check` green before allowing the merge, merge, tag (`vX.Y.Z`) on `main`. Update workflow triggers (`R-CMD-check.yaml`, `format-check.yaml`) to run on `dev` pushes and on PRs into `dev`/`main`. Update `.claude/CLAUDE.md` git section ("commit directly to main" -> "commit directly to dev; main only advances via release PR").
- **pkgdown site: GitHub Pages deployment**: Ship the pkgdown site with 1.0. `_pkgdown.yml` and local `pkgdown::build_site()` are already set up (vignettes written, README slimmed with links to them). Deployment is blocked on the repo being private -- once it's public (or `gh`/Pages access for private repos is sorted), run `usethis::use_pkgdown_github_pages()` to add the deploy workflow and enable Pages.

- check new iscn
- **[P2] Auto-fix missing comma after sex complement**: e.g. `47,XY+13[19]` (should be `47,XY,+13[19]`) glues the sex complement directly to the first abnormality with no separating comma. `validate_karyotypes()` splits clone 1 on commas, so the token becomes `XY+13`, matches neither `.sex_complements` nor the aberrant-sex-chromosome regex, and trips `no_sex_complement` -- currently unfixable since that issue type has no registered fix. Add a fix rule (like the `.dirty_patterns` entries in `R/assess.R`) that inserts a comma when a sex complement is immediately followed by `+`/`-`.
- **Fullwidth ISCN punctuation normalization**: `new_iscn.csv` has fullwidth brackets (`［`/`］`) used as literal metaphase-count bracket notation, plus other fullwidth punctuation (tilde, parens, comma, semicolon, equals). Currently lands in `stray_non_ascii` as unfixable. Needs its own normalization pattern in `R/assess.R` -- larger scope than the X/Y homoglyph fix, separate task.
- **`trailing_narrative` detect is Latin-script-biased**: its `detect` regexes require `[A-Z][a-z]`-shaped narrative, so CJK-only trailing narrative (seen in `new_iscn.csv`) is silently stripped by the fix but never reported by `check_karyo()` (confirmed via `.run_dirty_patterns()` in `R/assess.R`). Not a parsing bug -- the correct karyotype is still extracted -- but a reporting-accuracy gap worth a pass later.
- check that cols are propagated with every command (and even when id_col is not declared)
- make commands pipe friendly
- add a vignette that explain errors from check/preprocess with examples and potential solutions
- manually review tests with non-ASCII and other weird placeholders