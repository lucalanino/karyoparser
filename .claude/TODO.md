# TODO

## Before 1.0

### Parsing & data quality

- **[P3] Manually review tests with non-ASCII and other weird placeholders**.

### Release process

- **[P3] Branch strategy (decided, not yet implemented)**: `main` stays the default branch (so untagged `pak::pak()`/`install_github()` installs always resolve to the last release -- confirmed via pak docs: "if `<detail>` is missing, the latest commit of the default branch is used"). Ongoing work moves to a `dev` branch (not set as default). Releases: PR `dev -> main`, require `R-CMD-check` green before allowing the merge, merge, tag (`vX.Y.Z`) on `main`. Update workflow triggers (`R-CMD-check.yaml`, `format-check.yaml`) to run on `dev` pushes and on PRs into `dev`/`main`. Update `.claude/CLAUDE.md` git section ("commit directly to main" -> "commit directly to dev; main only advances via release PR").
- **[P3] pkgdown site: GitHub Pages deployment**: Ship the pkgdown site with 1.0. `_pkgdown.yml` and local `pkgdown::build_site()` are already set up (vignettes written, README slimmed with links to them). Deployment is blocked on the repo being private -- once it's public (or `gh`/Pages access for private repos is sorted), run `usethis::use_pkgdown_github_pages()` to add the deploy workflow and enable Pages.
