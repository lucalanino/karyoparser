# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Working with Claude

- After edits are made and explained, ask before committing/pushing — don't do it automatically. Skip this only if explicitly asked.
- After larger/multi-file edits, run `/simplify` to check for reuse, quality, and efficiency issues. Skip for trivial or doc-only changes.

## Project Overview

R package (`karyoparser`) for parsing ISCN karyotype strings into structured binary features. Focused on myeloid neoplasm-related chromosomal aberrations. Installable via `devtools::install_github("lucalanino/karyo-parser")`.

## Versioning

Scheme: `major.minor.patch` (semantic versioning). Current pre-release series: `0.x.y`. First production release will be `1.0.0`.

| Component | When to bump | Examples |
|---|---|---|
| `patch` | Bug fixes, internal refactors, test additions, doc-only changes | `0.1.0` → `0.1.1` |
| `minor` | New exported function, new aberration rule, new parameter, behavior change | `0.1.0` → `0.2.0` |
| `major` | Breaking API change after 1.0.0 (renamed/removed parameter or column, changed return structure). In the `0.x` series, breaking changes may be released as minor bumps. | `0.x.y` → `1.0.0` (first production release) |

No bump needed for: formatting-only commits, CI/tooling changes, README edits.

Three places to update:
1. `Version:` field in `DESCRIPTION`
2. `**Version**:` line in `README.md`
3. Version assertion in `tests/testthat/test_parse_karyo.R`

The `.karyoparser_version` constant in `R/parse_karyo-package.R` is derived from `DESCRIPTION` automatically — no manual update needed there.

## After any user-facing change

Run in order:
1. Review the roxygen `#'` block for any changed exported function — check `@param` defaults, `@return` column list, and caller attributions in any related internal-function docs.
2. `devtools::document()` — regenerate man/ pages
3. `devtools::check()` — confirm 0 errors/warnings/notes
4. Bump version in `DESCRIPTION` and `README.md`

## First-time setup (per machine)

```bash
# Point git to the tracked hooks directory
git config core.hooksPath .githooks
```

This enables the pre-push hook that runs `air format . --check` before every push, blocking unformatted code from reaching the remote.

## Development Commands

```bash
# Run tests
Rscript -e 'devtools::test()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test_parse_karyo.R")'

# Regenerate NAMESPACE and man/ pages after changing @export or roxygen docs
Rscript -e 'devtools::document()'

# Full R CMD check (0 errors/warnings/notes expected)
Rscript -e 'devtools::check()'

# Load package in dev mode (without installing)
Rscript -e 'devtools::load_all()'

# Format R source files with air
air format .
```

## Workflow design principle

The two stages of the workflow are **strictly separate** and must remain so:

1. **`preprocess_karyo()`** — the one and only place where cleaning and normalization happen. This includes dirty-pattern fixes, unicode normalization, whitespace collapsing, and delimiter tightening. Output is fully normalized and parser-ready.
2. **`parse_karyo()`** — pure parser. Trusts its input completely. Does **no** cleaning or normalization, not even unicode or whitespace. If called with `on_issues="fix"`, it routes input through `preprocess_karyo()` first (via `.assess_karyotypes()`). If called with `on_issues="warn"` or `"stop"`, it flags issues and warns/errors — it does not fix them.

**Never add cleaning or normalization logic to `parse_karyo()`.** Any new normalization step belongs in `preprocess_karyo()` and, if it represents a detectable data quality issue, should be added as a named entry in `.dirty_patterns` so `check_karyo()` can flag it.

## Architecture

Package source is split across multiple files in `R/`. A legacy standalone `source()`-able script is in `legacy/parse_karyo.R` (build-ignored).

| File | Contents |
|---|---|
| `R/parse_karyo-package.R` | `"_PACKAGE"`, `globalVariables()`, `.karyoparser_version`, `.sex_complements` |
| `R/rules.R` | `rules_table()` |
| `R/preprocess.R` | `.dirty_patterns`, `.sex_alt`, `.fixable_issue_types`, `.all_issue_types` (constants), `preprocess_karyo()`, `flag_unpreprocessed()` (internal), `normalize_iscn()` (internal) |
| `R/validate.R` | `.assess_karyotypes()` (internal), `.apply_dirty_fixes()` (internal), `check_karyo()`, `validate_karyotypes()` (internal) |
| `R/ploidy.R` | `ploidy_from_count()` (internal), `ploidy_category()` (internal), `extract_clone_data()` (internal) |
| `R/helpers.R` | `truncate_str()`, `empty_issues_tibble()`, `strip_bands()`, `normalize_token()`, `.karyotype_col_candidates`, `.id_col_candidates` (constants), `.extract_df_input()` (internal) |
| `R/parse_karyo.R` | `parse_karyo()` + internal pipeline helpers |

### Exported API

- **`parse_karyo()`** — Main entry point. Accepts character vector, data frame, or `karyo_preprocessed` tibble. When a `karyo_preprocessed` object is passed, `$preprocessed` is used automatically and cached indices + id column are reused — no extra arguments needed. Returns tibble with binary aberration flags, ploidy, monosomy/trisomy columns, derived flags (complex, monosomal, mixed ploidy), and `fixable_error`/`unfixable_error` (0L/1L) integer columns. `verbose = FALSE` default. Uses `on_issues` (`"fix"` / `"warn"` / `"stop"`). Default `"fix"` auto-applies `preprocess_karyo()` to dirty rows and truncates chimeric rows to the first clone before `//`. `"warn"` returns NA for all issue rows. `"stop"` errors immediately. Structural errors are always NA regardless.
- **`check_karyo(x, verbose = FALSE, karyotype_column = NULL, id_column = NULL)`** — Unified pre-parse check. Accepts character vector or data frame. For data frame input, karyotype column is auto-detected from `.karyotype_col_candidates`; id column is auto-detected from `.id_col_candidates` and placed first in the output. Returns a `karyo_check` tibble: one row per input, one 0L/1L column per issue type, plus `fixable` and `unfixable` summary columns. `verbose=FALSE` is total silence; `verbose=TRUE` prints count + per-issue breakdown and returns visibly. The result carries `.kp_assessment`, `.kp_raw`, `.kp_karyotype_col`, `.kp_id_col`, `.kp_id_values` as attributes.
- **`rules_table()`** — Regex rules defining which aberrations to detect. Priority system (100 > 95 > 90 > 85 > 60-80) resolves overlapping matches.
- **`preprocess_karyo(x, karyotype_column = NULL, id_column = NULL, verbose = FALSE)`** — The complete cleaning and normalization pipeline. Accepts a character vector, data frame, or `karyo_check` tibble. For data frame input: same column detection as `check_karyo()`. When given a `karyo_check`, reuses cached assessment and propagates id column. Fixes dirty markers then runs `normalize_iscn()`. Returns a `karyo_preprocessed` tibble with optional id column (first), `original`, `preprocessed`, and `status` columns. Carries `.kp_fixable_rows`, `.kp_unfixable_rows`, `.kp_chimeric_rows`, `.kp_chimeric_all_rows`, `.kp_id_col`, `.kp_id_values` as attributes. Passing it directly to `parse_karyo()` requires no extra arguments.

### Internal pipeline helpers (in `R/parse_karyo.R`, unexported)

- `build_sample_meta()` — Bracket parsing, ploidy classification, normal karyotype detection.
- `build_clone_tokens()` — Clone splitting, token expansion, normalization.
- `match_rules()` — Regex matching via pre-allocated logical matrix with priority-based deduplication.
- `compute_aneuploidy()` — Monosomy/trisomy flag detection.
- `compute_comma_counts()` — Idem-aware comma-based aberration counting.
- `compute_unique_counts()` — Unique aberration counting (for complex flag).
- `blank_rows()` — Vectorized NA-row generation for invalid inputs.

### Other internal helpers (unexported)

- `truncate_str()`, `empty_issues_tibble()`, `strip_bands()`, `normalize_token()`, `.extract_df_input()` (in `R/helpers.R`)
- `.karyotype_col_candidates`, `.id_col_candidates` constants (in `R/helpers.R`)
- `flag_unpreprocessed()`, `normalize_iscn()` (in `R/preprocess.R`)
- `.apply_dirty_fixes()`, `validate_karyotypes()` (in `R/validate.R`)
- `ploidy_from_count()`, `ploidy_category()`, `extract_clone_data()` (in `R/ploidy.R`)

### Key design details

- **`.pk_row_id`**: Internal row-tracking column, deliberately named to avoid collision with user data columns like `sample_id`. Dropped before output.
- **`.sex_complements`**: Package-level constant listing valid sex chromosome complements. Used in both validation and parsing.
- **`.dirty_patterns`**: Named list of pre-normalization issue patterns in `R/preprocess.R`. Each entry has `detect`, `use_trimmed`, `fix` (empty list = detectable but not fixable), and `detail`. Entries: `unicode_notation`, `embedded_newline`, `html_entities`, `leading_dot`, `fish_notation`, `trailing_narrative`, `midstring_linewrap`, `missing_sex_comma`, `mar_space`, `zero_host_chimera`. After all dirty fixes are applied, `normalize_iscn()` runs as a final normalization pass (whitespace collapsing, delimiter tightening, `\s*\[` tightening).
- **`.fixable_issue_types`**: Constant in `R/preprocess.R`. Dirty pattern names with non-empty `fix` lists + `"chimeric_separator"`.
- **`.all_issue_types`**: Constant in `R/preprocess.R`. Fixed schema for `check_karyo()` columns — all detectable issue types in display order.
- **`.apply_dirty_fixes()`**: Private function in `R/validate.R`. Single source of truth for the dirty-fix loop (applies all `.dirty_patterns` fixes). Called by `.assess_karyotypes()` (Step 2).
- **Priority-based rule matching**: A logical matrix (tokens × rules) is built with one `str_detect` call per rule. For each token, the highest-priority matching rule is selected via `which.max`. Ties across tokens sharing the same `aberr_norm` within a clone are resolved by `slice_head` after sorting by priority desc.
- **Idem expansion**: For `comma_count_aberrations`, clones with "idem" inherit the stemline's (clone 1) aberration count. Both "idem" and "sl" are excluded from unique aberration counts used for `complex_karyotype`.
- **Monosomal karyotype**: Requires ≥2 autosomal monosomies OR ≥1 autosomal monosomy + ≥1 structural aberration where `counts_for_monosomal=TRUE`. Sex chromosome monosomies and `marker_chromosome` don't qualify.
- **Ploidy gaps**: 34-39 and 47-50 map to "other" intentionally (not classified as hypo/hyperdiploid).
- **Composite karyotypes**: For ploidy, only clones with ≥5 metaphases are eligible (when brackets exist). Most abnormal clone (furthest from 46) determines classification.
- **Deduplication**: When input contains duplicate karyotype strings, only unique strings are parsed; results are joined back to all rows. Zero overhead when no duplicates exist.

### `parse_karyo()` pipeline flow

```
Input → [Input routing: detect is_preprocessed FIRST]:
    karyo_preprocessed: raw_vec = $preprocessed, original_vec = $original, id from attrs
    plain data frame:   .extract_df_input() → raw_vec, original_vec, id
    character vector:   raw_vec = original_vec = as.character(karyotypes)
  [Guard]:
    is_preprocessed: reuse cached .kp_* indices, skip .assess_karyotypes()
    standard path:   .assess_karyotypes(raw_vec) → indices derived
    "stop": error if any issues; raw_vec = normalize_iscn(raw_vec)
    "fix":  raw_vec = $processed (NA for unfixable); unfixable → NA
    "warn": issue rows → NA (no fixing); raw_vec = normalize_iscn(raw_vec)
  Build input_df (original_karyotype = original_vec, normalized_karyotype = raw_vec)
  Filter issue rows → [Dedup] → Parse unique:
    build_sample_meta → build_clone_tokens →
    match_rules / compute_aneuploidy / compute_comma_counts / compute_unique_counts →
    Assembly (monosomal/complex inline) →
  [Join back if deduped] → Recombine issue rows → Attach fixable_error/unfixable_error →
  Attach id column (from id_values/id_col_name) → Output
```

Note: `normalized_karyotype` column = `normalize_iscn(raw)` in all three modes. In `"fix"` mode, `raw_vec` is also dirty-fixed and chimeric-truncated before normalization. `fixable_error = 1` for any row that had a fixable issue (dirty marker or chimeric separator), regardless of whether it was fixed. When input is `karyo_preprocessed`, `original_karyotype` = `$original` (the truly raw string), not the preprocessed form.

`normalize_iscn()` is called inside `preprocess_karyo()` and (for "warn"/"stop") inside `parse_karyo()`. `parse_karyo()` applies no other normalization.

### Adding new aberration rules

Add rows to `rules_table()` in `R/rules.R`: `flag_name` (output column), `regex`, `category`, `priority`, `counts_for_monosomal`. Run `devtools::document()` after.

## Testing

539 assertions (sections 1–29) in `tests/testthat/test_parse_karyo.R` covering: all regex rules (positive/negative/reversed), priority system, ploidy classification (incl. gap boundaries, ineligible clones, range averaging), monosomy/trisomy detection, complex/monosomal flags, CBF-AML override (t(8;21)/inv(16)/t(16;16)), preprocessing, idem expansion (incl. clean-stemline edge case), `check_karyo()`, `on_issues` guard (`"fix"`, `"warn"`, `"stop"`), `fixable_error`/`unfixable_error` columns, ID column detection (incl. explicit `karyotype_column`), deduplication, multi-group rule firing, `preprocess_karyo()`, `.dirty_patterns`, trailing narrative (4-rule chain), midstring_linewrap (incl. `+`), fish_notation, mar_space, missing_sex_comma, chimeric_separator, updated_iscn, zero_host_chimera, unicode_notation, embedded_newline, normalize_iscn (via preprocess path), `rules_table()` structure, `normalized_karyotype` column, workflow consistency (check → preprocess → parse chain, false-positive prevention, fixable_error cross-mode consistency), entry-point consistency (data frame input for all three functions, id propagation through full pipeline, type guard, `karyo_preprocessed` without `karyotype_column`, verbose alignment), and edge cases.

Tests use a `pk()` helper that wraps `parse_karyo(..., on_issues = "warn", verbose = FALSE)`.

## TODO — Before 1.0 Release

- **API audit**: Review every exported name before cutting 1.0 — after that, renames are breaking changes. Cover:
  - Exported function names: `parse_karyo()`, `check_karyo()`, `preprocess_karyo()`, `rules_table()`
  - Parameters: `on_issues = c("fix","warn","stop")` — names and defaults
  - Output columns of `parse_karyo()`: all binary aberration flags, `fixable_error`/`unfixable_error`, `original_karyotype`, `normalized_karyotype`, `ploidy`, etc.
  - Output columns of `check_karyo()`: column order, `fixable`/`unfixable` naming

## TODO - Future Discussions

- **pkgdown site + GitHub Pages**: Set up a documentation website via `usethis::use_pkgdown()`, hosted on GitHub Pages. Plan:
  - Convert `README.md` to `README.Rmd` — source file with executable R code chunks; `devtools::build_readme()` regenerates `README.md` from it. Use `dplyr::select()` or `knitr::kable()` to show curated, narrow output (not the full wide tibble with all flag columns). Never hand-edit `README.md` once `README.Rmd` exists.
  - Trim the README to a lean intro (what it does, install snippet, one good example) — it becomes the site homepage
  - Write vignettes in `vignettes/` for full docs: get-started, aberration rules, preprocessing/QC, ploidy
  - Function reference pages are auto-generated from roxygen docs
  - **CRITICAL**: Audit `.Rbuildignore` and `.gitignore` before building — `pkgdown::build_site()` can inadvertently pull gitignored files (e.g. `dev-data/`) into `docs/`. Verify nothing sensitive leaks into the built site before pushing.
- **CI workflow (GitHub Actions)**: `usethis::use_github_actions("check-standard")` runs `R CMD check` automatically on every push across platforms. Also useful as a guard against accidentally committing ignored files or broken states. Note: free tier is 2,000 min/month on private repos (~1,000 pushes at ~2 min/check) — monitor if it becomes a constraint. **Add coverage upload to the same workflow** via `covr::codecov()` + Codecov GitHub App — do both together, not separately.
- **Spell check**: `usethis::use_spell_check()` adds spell checking of roxygen docs and README to the test suite.
- **Discuss: per-clone export option**: Option to split composite karyotypes into one row per clone in the output, instead of one row per ISCN string. Discuss API design, how derived flags (complex, monosomal) behave per-clone vs. per-karyotype, and whether this is a parameter on `parse_karyo()` or a post-processing helper. Each clone row should include a `clone_abundance` column (percentage of total metaphases belonging to that clone, derived from bracket counts).
