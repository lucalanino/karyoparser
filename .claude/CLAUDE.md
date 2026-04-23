## Project Overview

`karyoparser` is an R package for parsing ISCN karyotype strings into structured features for downstream analysis.

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

### Testing

- Tests for `R/{name}.R` go in `tests/testthat/test-{name}.R`. 
- All new code should have an accompanying test.
- If there are existing tests, place new tests next to similar existing tests.
- Strive to keep your tests minimal with few comments.

### Documentation

- Every user-facing function should be exported and have roxygen2 documentation at the top of the script.
- Wrap roxygen comments at 80 characters.
- Internal functions should not have roxygen documentation.
- Always re-document the package after changing a roxygen2 comment.

### Versioning

- Follow the scheme: `major.minor.patch`
- No bump needed for: formatting-only commits, CI/tooling changes, README-only edits.

### TODOs

- Tracked in `.claude/TODO.md`. Entries should be brief, informative, and prioritized.

### Core Architecture

- Surfaced functions should have one job: checking, preprocessing or parsing