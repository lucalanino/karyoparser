# TODO

## Before 1.0

### Parsing & data quality

- **[P2] Myeloid rule completeness**: Final pass over `myeloid_rules`.
- **[P2] Auto-fix missing comma after sex complement glued to a letter-token**: e.g. `46,XYdel(5q)` (should be `46,XY,del(5q)`) -- the zero-separator counterpart of the `+`/`-` case already fixed in `missing_sex_comma`. The pre-existing space-glued variant (`46,XY del(5q)`) already works; only the no-space form is unhandled, and still trips `no_sex_complement` as unfixable. Riskier than the `+`/`-` case: `[a-z(]` right after a sex complement is a much broader trigger surface, so needs real-corpus validation (same process as the `+`/`-` fix) to check for false positives before adding a detect/fix pair to `missing_sex_comma` in `R/assess.R`, plus a check that it can't collide with the constitutional-sex-complement `c` suffix (e.g. `47,XXYc`).
- **[P2] Fullwidth ISCN punctuation normalization**: some real-world karyotype strings use fullwidth brackets (`［`/`］`) as literal metaphase-count bracket notation, plus other fullwidth punctuation (tilde, parens, comma, semicolon, equals). Currently lands in `stray_non_ascii` as unfixable. Needs its own normalization pattern in `R/assess.R` -- larger scope than the X/Y homoglyph fix, separate task. Implement the fix as a single `stringr::str_replace_all(x, c(pattern = replacement, ...))` call with a named vector rather than one `list(pattern=, replacement=)` entry per character -- benchmarked ~2x faster (one `stri_replace_all_regex(..., vectorize_all = FALSE)` call vs. N sequential `str_replace_all()` calls) and more compact.
- **[P2] Finish real-world-data review**: continue auditing real-world karyotype strings for unhandled dirty-data patterns, beyond the fullwidth-punctuation finding above (missing sex-complement comma was fixed separately) and the `trailing_narrative` bias below.
- **[P3] `trailing_narrative` detect is Latin-script-biased**: its `detect` regexes require `[A-Z][a-z]`-shaped narrative, so CJK-only trailing narrative is silently stripped by the fix but never reported by `check_karyo()` (confirmed via `.run_dirty_patterns()` in `R/assess.R`). Not a parsing bug -- the correct karyotype is still extracted -- but a reporting-accuracy gap worth a pass later.
- **[P3] Manually review tests with non-ASCII and other weird placeholders**.

### API ergonomics

- **[P2] Column propagation audit**: check that columns are propagated correctly by every command, including when `id_col` is not declared.
- **[P2] Pipe-friendly commands**: make commands pipe friendly.
- **[P3] Reconsider `check_karyo()`/`preprocess_karyo()` output shape**: the fixable-issue-type list keeps growing (unicode, dots, separators, glued commas, ...) and auto-fixing is becoming a first-class feature rather than an edge case. Current output is one wide 0/1 column per issue type with no visibility into *what changed*. Worth exploring something that shows before/after evidence per row (e.g. a diff or fixed-span annotation), not just a flag.

### Performance

- **[P3] Speed up `.dirty_patterns` fixing in `R/assess.R`**: `.run_dirty_patterns()` applies each pattern's `fix` list as a loop of sequential `stringr::str_replace_all()` calls (one per `list(pattern=, replacement=)` entry); collapsing same-pass, non-overlapping char/token substitutions (e.g. `unicode_notation`'s dash/space variants) into a single named-vector `str_replace_all()` call benchmarked ~2x faster (see `fullwidth_punctuation` above). Worth a pass over the other multi-entry `fix` lists to see which ones qualify, and whether row-batching elsewhere in the check/preprocess/parse pipeline has similar wins.
- **[P3] Profile each stage of `check_karyo()`/`preprocess_karyo()`/`parse_karyo()`**: benchmark the major internal steps (dirty-pattern fixing, validation, tokenization, translocation/aberration flagging, output-schema assembly, etc.) on a realistic-size input to see where time actually goes, before sinking effort into speedups. Prioritize by measured share of runtime, not by guesswork -- some sections may already be fast enough that optimizing them is a waste of time.

### Release process

- **[P3] Branch strategy (decided, not yet implemented)**: `main` stays the default branch (so untagged `pak::pak()`/`install_github()` installs always resolve to the last release -- confirmed via pak docs: "if `<detail>` is missing, the latest commit of the default branch is used"). Ongoing work moves to a `dev` branch (not set as default). Releases: PR `dev -> main`, require `R-CMD-check` green before allowing the merge, merge, tag (`vX.Y.Z`) on `main`. Update workflow triggers (`R-CMD-check.yaml`, `format-check.yaml`) to run on `dev` pushes and on PRs into `dev`/`main`. Update `.claude/CLAUDE.md` git section ("commit directly to main" -> "commit directly to dev; main only advances via release PR").
- **[P3] pkgdown site: GitHub Pages deployment**: Ship the pkgdown site with 1.0. `_pkgdown.yml` and local `pkgdown::build_site()` are already set up (vignettes written, README slimmed with links to them). Deployment is blocked on the repo being private -- once it's public (or `gh`/Pages access for private repos is sorted), run `usethis::use_pkgdown_github_pages()` to add the deploy workflow and enable Pages.
