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
* No non-ASCII characters in R source files (CRAN requirement). In comments use `--` and `->` instead of `—`/`→`. In string literals use `\uXXXX` escapes (e.g. `—` for em dash).
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

- Unless explicitly told otherwise, commit directly to `main`.

### TODOs

- Tracked in `.claude/TODO.md`. Entries should be brief, informative, and prioritized.
