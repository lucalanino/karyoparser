# TODO

## Before 1.0

### Parsing & data quality

- **[P2] Myeloid rule completeness**: Final pass over `myeloid_rules`.
- **[P2] Auto-fix missing comma after sex complement**: e.g. `47,XY+13[19]` (should be `47,XY,+13[19]`) glues the sex complement directly to the first abnormality with no separating comma. `validate_karyotypes()` splits clone 1 on commas, so the token becomes `XY+13`, matches neither `.sex_complements` nor the aberrant-sex-chromosome regex, and trips `no_sex_complement` -- currently unfixable since that issue type has no registered fix. Add a fix rule (like the `.dirty_patterns` entries in `R/assess.R`) that inserts a comma when a sex complement is immediately followed by `+`/`-`.
- **[P2] Fullwidth ISCN punctuation normalization**: `new_iscn.csv` has fullwidth brackets (`［`/`］`) used as literal metaphase-count bracket notation, plus other fullwidth punctuation (tilde, parens, comma, semicolon, equals). Currently lands in `stray_non_ascii` as unfixable. Needs its own normalization pattern in `R/assess.R` -- larger scope than the X/Y homoglyph fix, separate task. Implement the fix as a single `stringr::str_replace_all(x, c(pattern = replacement, ...))` call with a named vector rather than one `list(pattern=, replacement=)` entry per character -- benchmarked ~2x faster (one `stri_replace_all_regex(..., vectorize_all = FALSE)` call vs. N sequential `str_replace_all()` calls) and more compact.
- **[P2] Finish `new_iscn.csv` review**: Broader real-world-data pass over `new_iscn.csv` that already produced the two findings above (missing sex-complement comma, fullwidth punctuation) plus the `trailing_narrative` bias below; keep scanning the rest of the file for other unhandled patterns.
- **[P3] `trailing_narrative` detect is Latin-script-biased**: its `detect` regexes require `[A-Z][a-z]`-shaped narrative, so CJK-only trailing narrative (seen in `new_iscn.csv`) is silently stripped by the fix but never reported by `check_karyo()` (confirmed via `.run_dirty_patterns()` in `R/assess.R`). Not a parsing bug -- the correct karyotype is still extracted -- but a reporting-accuracy gap worth a pass later.
- **[P3] Manually review tests with non-ASCII and other weird placeholders**.

### API ergonomics

- **[P2] Column propagation audit**: check that columns are propagated correctly by every command, including when `id_col` is not declared.
- **[P2] Pipe-friendly commands**: make commands pipe friendly.

### Performance

- **[P3] Speed up `.dirty_patterns` fixing in `R/assess.R`**: `.run_dirty_patterns()` applies each pattern's `fix` list as a loop of sequential `stringr::str_replace_all()` calls (one per `list(pattern=, replacement=)` entry); collapsing same-pass, non-overlapping char/token substitutions (e.g. `unicode_notation`'s dash/space variants) into a single named-vector `str_replace_all()` call benchmarked ~2x faster (see `fullwidth_punctuation` above). Worth a pass over the other multi-entry `fix` lists to see which ones qualify, and whether row-batching elsewhere in the check/preprocess/parse pipeline has similar wins.
- **[P3] Profile each stage of `check_karyo()`/`preprocess_karyo()`/`parse_karyo()`**: benchmark the major internal steps (dirty-pattern fixing, validation, tokenization, translocation/aberration flagging, output-schema assembly, etc.) on a realistic-size input to see where time actually goes, before sinking effort into speedups. Prioritize by measured share of runtime, not by guesswork -- some sections may already be fast enough that optimizing them is a waste of time.

### Documentation

- **[P3] Errors vignette**: add a vignette explaining errors from `check_karyo()`/`preprocess_karyo()` with examples and potential solutions.

### Release process

- **[P3] Branch strategy (decided, not yet implemented)**: `main` stays the default branch (so untagged `pak::pak()`/`install_github()` installs always resolve to the last release -- confirmed via pak docs: "if `<detail>` is missing, the latest commit of the default branch is used"). Ongoing work moves to a `dev` branch (not set as default). Releases: PR `dev -> main`, require `R-CMD-check` green before allowing the merge, merge, tag (`vX.Y.Z`) on `main`. Update workflow triggers (`R-CMD-check.yaml`, `format-check.yaml`) to run on `dev` pushes and on PRs into `dev`/`main`. Update `.claude/CLAUDE.md` git section ("commit directly to main" -> "commit directly to dev; main only advances via release PR").
- **[P3] pkgdown site: GitHub Pages deployment**: Ship the pkgdown site with 1.0. `_pkgdown.yml` and local `pkgdown::build_site()` are already set up (vignettes written, README slimmed with links to them). Deployment is blocked on the repo being private -- once it's public (or `gh`/Pages access for private repos is sorted), run `usethis::use_pkgdown_github_pages()` to add the deploy workflow and enable Pages.
