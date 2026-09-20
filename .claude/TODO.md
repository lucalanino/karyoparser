# TODO

## After 1.0

### Parsing & data quality


### Release process

- **[P2] Branch strategy (decided, partly implemented)**: `main` stays the default branch, so untagged `pak::pak()`/`install_github()` installs resolve to the last release -- confirmed via pak docs: "if `<detail>` is missing, the latest commit of the default branch is used". `.claude/CLAUDE.md` git section is already updated. Remaining: create `dev` from `main` at the next change and work there; releases go PR `dev -> main`, require `R-CMD-check` green, merge, then annotated tag `vX.Y.Z` on `main`. Update workflow triggers -- `R-CMD-check.yaml` (add `dev` to both `push` and `pull_request`) and `format-check.yaml` (add `dev` to `push`, and add a `pull_request` trigger, which it currently lacks entirely, so release PRs into `main` get no format check). `pkgdown.yaml` likely needs no change, since its stock `push: branches: [main]` trigger then fires only on release merges -- confirm rather than assume. Optional later, deliberately skipped for now: branch protection on `main`, `NEWS.md`, and the `x.y.z.9000` dev-version convention.
