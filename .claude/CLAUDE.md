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
- No bump needed for: formatting-only commits, CI/tooling changes, README-only edits.

### Git

- Always check the current branch before committing.
- All changes live on `dev`. If `dev` does not exist yet, create it from
  `main` (`git switch -c dev`) rather than committing to `main`.
- `main` only advances via a PR from `dev`, merged once `R-CMD-check` is
  green. Never commit directly to `main` unless explicitly told otherwise.
- Releases are annotated tags (`vX.Y.Z`) on `main`.

### TODOs

- Tracked in `.claude/TODO.md`. Entries should be brief, informative, and prioritized.
