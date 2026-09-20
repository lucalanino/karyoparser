# TODO

## Before 1.0

### Parsing & data quality

- **[P2] Myeloid rule completeness**: Final pass over `myeloid_rules`.
- **[P2] Finish real-world-data review**: continue auditing real-world karyotype strings for unhandled dirty-data patterns (missing sex-complement comma and fullwidth ISCN punctuation were both fixed separately) and the `trailing_narrative` bias below.
- **[P3] Manually review tests with non-ASCII and other weird placeholders**.

### Risk stratification

- **[P2] ELN 2022 cytogenetic risk (`assign_risk(scheme = "eln2022")`)**: IPSS-R
  shipped first; ELN follows on the same API. Specifics to get right:
  - **Complex uses the ELN carve-out, not the generic `complex_karyotype`**:
    excludes hyperdiploid karyotypes with 3+ trisomies/polysomies and no
    structural abnormality. `complex_karyotype` stays scheme-neutral because
    IPSS-R has no such exclusion.
  - **KMT2A carve-out**: `t(9;11)(p21;q23)` is Intermediate and takes
    precedence over the rarer adverse lesions; all *other* KMT2A/11q23
    rearrangements (`t_v_11q23`) are Adverse. Order these so t(9;11) wins.
  - **Class-defining favorable lesions override adverse-complex**, mirroring
    the CBF-AML override already applied to `monosomal_karyotype`.
  - Name it `eln2022_cyto_risk`: ELN 2022 is a *genetic* classification
    needing NPM1/FLT3-ITD/CEBPA/TP53/ASXL1/RUNX1 etc. Cytogenetics alone
    cannot assign a true ELN risk group.
- **[P3] Risk vignette**: fold both schemes into one article once ELN lands,
  leading with the "cytogenetic component only" caveat.

### Release process

- **[P3] Branch strategy (decided, not yet implemented)**: `main` stays the default branch (so untagged `pak::pak()`/`install_github()` installs always resolve to the last release -- confirmed via pak docs: "if `<detail>` is missing, the latest commit of the default branch is used"). Ongoing work moves to a `dev` branch (not set as default). Releases: PR `dev -> main`, require `R-CMD-check` green before allowing the merge, merge, tag (`vX.Y.Z`) on `main`. Update workflow triggers (`R-CMD-check.yaml`, `format-check.yaml`) to run on `dev` pushes and on PRs into `dev`/`main`. Update `.claude/CLAUDE.md` git section ("commit directly to main" -> "commit directly to dev; main only advances via release PR").
- **[P3] pkgdown site: GitHub Pages deployment**: Ship the pkgdown site with 1.0. `_pkgdown.yml` and local `pkgdown::build_site()` are already set up (vignettes written, README slimmed with links to them). Deployment is blocked on the repo being private -- once it's public (or `gh`/Pages access for private repos is sorted), run `usethis::use_pkgdown_github_pages()` to add the deploy workflow and enable Pages.
