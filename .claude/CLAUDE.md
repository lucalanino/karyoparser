## Project Overview

`karyoparser` is an R package for parsing ISCN karyotype strings into tabular format.


### Key commands

```
# To run code
Rscript -e "devtools::load_all(); code"

# To run all tests
Rscript -e "devtools::test()"

# To run all tests for files starting with {name}
Rscript -e "devtools::test(filter = '^{name}')"

# To run all tests for R/{name}.R
Rscript -e "devtools::test_active_file('R/{name}.R')"

# To run a single test "blah" for R/{name}.R
Rscript -e "devtools::test_active_file('R/{name}.R', desc = 'blah')"

# To redocument the package
Rscript -e "devtools::document()"

# To reknit README.md from README.Rmd
Rscript -e "devtools::build_readme()"

# To check pkgdown documentation
Rscript -e "pkgdown::check_pkgdown()"

# To check the package with R CMD check
Rscript -e "devtools::check()"

# To format code
air format .
```

### Coding

* Always run `air format .` after generating code
* Use the base pipe operator (`|>`) not the magrittr pipe (`%>%`)
* Don't use `_$x` or `_$[["x"]]` since this package must work on R 4.1.
* Use `\() ...` for single-line anonymous functions. For all other cases, use `function() {...}`
* No non-ASCII characters in R source files (CRAN requirement). In comments and
  roxygen use `--` and `->` rather than the em dash and arrow glyphs.
* Never use an em dash as punctuation anywhere -- not in code, comments, roxygen,
  messages, tests, vignettes or README -- and not as a `\u2014` escape either.
  Use `--`. This is a style rule, not a portability one: the escape form is
  perfectly CRAN-safe, it is just not wanted.
* `\uXXXX` escapes remain correct for non-ASCII characters the package must match
  in *data*. The dash character classes in `.dirty_patterns$unicode_notation`
  (`R/assess.R`) and `normalize_iscn()` (`R/preprocess_karyo.R`) contain `\u2014`
  on purpose: the package detects an em dash in dirty karyotype input and
  normalizes it to `-`. Do not strip those -- removing them silently disables
  the fix.
* Section header comments (rare -- most files are small and don't need them) use the RStudio/VS Code outline-navigable form `# Section Name ----`, one level (`# `) by default; nest `## Subsection Name ----` only when a section genuinely has sub-groupings. Never use full-width `# ---- Name ----------` banners.

### Testing

- Tests for `R/{name}.R` go in `tests/testthat/test-{name}.R`. 
- All new code should have an accompanying test.
- If there are existing tests, place new tests next to similar existing tests.
- Strive to keep your tests minimal with few comments.

### Documentation

- Every user-facing function should be exported and have roxygen2 documentation at the top of the script.
- Wrap roxygen comments at 80 characters.
- Internal functions will not have roxygen documentation.
- Always re-document the package after changing a roxygen2 comment.
- `README.md` is generated from `README.Rmd` with `devtools::build_readme()`.

### Versioning

- Follow the scheme: `major.minor.patch`
- Bump once per release cycle, not per commit. No `.9000` dev versions.
- No bump needed for: formatting-only commits, CI/tooling changes, README-only edits.
- `NEWS.md`: add a bullet for each user-visible change under the top
  `# karyoparser x.y.z` heading, in the same commit as the change. The heading
  must match `Version:` in `DESCRIPTION`. Releases before 1.0.1 have no entry.

### Git

- Always check the current branch before committing.
- `dev` is a throwaway branch, one per release cycle. All changes live on it,
  never directly on `main`.
- Start a cycle from an up-to-date `main`:
  `git switch main && git pull && git switch -c dev`.
- Release steps, once the work on `dev` is done:
  1. Bump `Version:` and check the top `NEWS.md` heading matches.
  2. `gh pr create --base main --fill`; wait for green checks
     (`gh pr checks --watch`).
  3. `gh pr merge --merge --delete-branch` (merge commit only, so the
     individual commits stay on `main`; squash and rebase are disabled).
  4. `git switch main && git pull && git branch -d dev`.
  5. `git tag -a vX.Y.Z -m "karyoparser X.Y.Z" && git push origin vX.Y.Z`.
  6. `gh release create vX.Y.Z --title "karyoparser X.Y.Z"` with that
     version's `NEWS.md` section as the notes.
- `main` is protected by a ruleset: PR required, `R-CMD-check` and
  `format-check` must pass, no force-push or deletion.
- Tags are annotated, on `main`, and match `DESCRIPTION`. Never move or reuse a
  pushed tag; release a new patch instead.

### TODOs

- Tracked in `.claude/TODO.md`. Entries should be brief, informative, and prioritized.
