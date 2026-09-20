# Assign Cytogenetic Risk Categories

Adds a cytogenetic risk category to each row of a
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
result. Interpretation only – no parsing or cleaning happens here; the
function reads columns
[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md)
already produced.

## Usage

``` r
assign_risk(parsed, scheme = c("ipssr", "eln2022"))
```

## Arguments

- parsed:

  A tibble returned by
  [`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md).
  Must carry the default
  [myeloid_rules](https://lucalanino.github.io/karyoparser/reference/myeloid_rules.md)
  flags; if `columns` was used to narrow the output, the classes the
  scheme needs must still be present.

- scheme:

  Character string naming the risk scheme: `"ipssr"` (default, IPSS-R
  cytogenetic risk group) or `"eln2022"` (ELN 2022 cytogenetic risk
  category). One per call.

## Value

`parsed` with the scheme's columns appended.

For `"ipssr"`:

- `ipssr_cyto_risk`: ordered factor, `Very Good` \< `Good` \<
  `Intermediate` \< `Poor` \< `Very Poor`; `NA` for rows that cannot be
  scored

- `ipssr_cyto_score`: integer 0-4 matching the category, `NA` alongside
  an `NA` category

For `"eln2022"`:

- `eln2022_cyto_risk`: ordered factor, `Favorable` \< `Intermediate` \<
  `Adverse`; `NA` for rows the parser could not score

## Details

**This is the cytogenetic component only, never a complete risk score.**
IPSS-R stratifies on five components – cytogenetics, marrow blast
percentage, haemoglobin, platelet count, and absolute neutrophil count –
of which this computes the first. ELN 2022 is a *genetic*
classification: assigning it properly needs *NPM1*, *FLT3*-ITD, *CEBPA*
bZIP, *TP53*, *ASXL1*, *RUNX1* and others, none of which a karyotype
carries. A normal karyotype is intermediate cytogenetically but
favorable with mutated *NPM1* and no *FLT3*-ITD, or adverse with mutated
*TP53*. Treat both columns as one input to a risk assessment, not as the
assessment.

Apply one scheme per call; chain for both:
`assign_risk(x, "ipssr") |> assign_risk("eln2022")`.

## Counting aberrations

"Isolated", "double", and "complex" are resolved with
`distinct_aberrations` – distinct normalized aberration tokens pooled
across all clones, which is the same count `complex_karyotype`
thresholds at 3. Rows the parser could not score (`distinct_aberrations`
is `NA`) return `NA`.

## IPSS-R categories

Following the IPSS-R cytogenetic risk groups (Schanz et al. 2012):

|  |  |  |
|----|----|----|
| Category | Score | Definition |
| Very Good | 0 | isolated `-Y`, isolated `del(11q)` |
| Good | 1 | normal, isolated `del(5q)`/`del(12p)`/`del(20q)`, double including `del(5q)` |
| Intermediate | 2 | `del(7q)`, `+8`, `+19`, `i(17q)`, any other single or double |
| Poor | 3 | `-7`, `inv(3)`/`t(3q)`/`del(3q)`, double including `-7`/`del(7q)`, 3 abnormalities |
| Very Poor | 4 | more than 3 abnormalities |

Two resolutions are applied where the published table is silent, both
deliberate and worth confirming against the source before relying on
them:

- A double containing **both** `del(5q)` and `-7`/`del(7q)` satisfies
  the Good row and the Poor row at once. Poor wins here, on the grounds
  that the more adverse lesion should not be masked by its partner.

- `t(3q)` is read as 3q26 rearrangement (`t_3q26_v`), not any 3q
  breakpoint, since the IPSS-R row reflects *MECOM*/EVI1 biology.
  `del(3q)` is matched at any 3q breakpoint.

A karyotype with no scoreable aberration token that is nonetheless not
normal – `45,X`, where the loss rides on the sex complement and `-X`
cannot be distinguished from `-Y` – returns `NA` rather than guessing.

## ELN 2022 categories

The cytogenetic rows of the ELN 2022 risk classification, evaluated in
strict precedence order:

1.  **Favorable** – `t(8;21)(q22;q22)`/*RUNX1::RUNX1T1*,
    `inv(16)(p13q22)` or `t(16;16)(p13;q22)`/*CBFB::MYH11*. These are
    class-defining, so they win over every adverse criterion below,
    mirroring the CBF-AML override already applied to
    `monosomal_karyotype`.

2.  **Intermediate, by carve-out** – `t(9;11)(p21;q23)`/*MLLT3::KMT2A*.
    Checked before the adverse rows because a `t(9;11)` also satisfies
    the adverse `t(v;11q23)`/*KMT2A*-rearranged row; ELN gives `t(9;11)`
    precedence, so it must be resolved first or every `t(9;11)` would be
    miscalled adverse.

3.  **Adverse** – `t(6;9)`/*DEK::NUP214*;
    `t(v;11q23)`/*KMT2A*-rearranged; `t(9;22)`/*BCR::ABL1*;
    `t(8;16)`/*KAT6A::CREBBP*; `inv(3)(q21q26)` or `t(3;3)(q21;q26)`;
    `t(3q26;v)`/*MECOM*-rearranged; `-5` or `del(5q)`; `-7`;
    `-17`/abn(17p); complex karyotype; monosomal karyotype.

4.  **Intermediate** – anything else, including a normal karyotype.

`t(6;9)` and `t(9;22)` are matched at both canonical and non-canonical
breakpoints, since ELN names the fusion rather than the breakpoint. This
matters: ELN writes `t(6;9)(p23.3;q34.1)`, which is a *non*-canonical
spelling for
[myeloid_rules](https://lucalanino.github.io/karyoparser/reference/myeloid_rules.md)
and would otherwise be missed. `inv(3)` and `t(3;3)` are matched only at
the canonical `q21q26`, because the non-canonical variants are not the
*GATA2*/*MECOM* lesion.

**Complex** uses the ELN definition, not the scheme-neutral
`complex_karyotype` column: three or more distinct aberrations,
excluding hyperdiploid karyotypes with three or more trisomies (or
polysomies) and no structural abnormality. The "in the absence of a
class-defining recurring genetic abnormality" clause falls out of the
precedence order – favorable lesions and `t(9;11)` are resolved before
complex is considered. `complex_karyotype` itself is deliberately left
alone, since IPSS-R's complex rows carry no such exclusion.

Two judgment calls, both worth confirming against the source:

- **abn(17p)** is read as `-17`, `del(17p)`, `add(17p)` or `i(17q)`.
  Including `i(17q)` reflects that it entails 17p loss.

- **Partial 17p loss implied by an unbalanced derivative is *not*
  counted**, because `unbal_partial_loss_*` is documented as separate
  from `del()`/`mono*` and folding it in here would breach that
  boundary. If you want those rows treated as abn(17p), combine
  `unbal_partial_loss_17p` yourself.

ELN 2022 has no numeric score, so `scheme = "eln2022"` adds only a
category column.

## See also

[`parse_karyo()`](https://lucalanino.github.io/karyoparser/reference/parse_karyo.md),
[`vignette("risk-stratification")`](https://lucalanino.github.io/karyoparser/articles/risk-stratification.md)

## Examples

``` r
parsed <- parse_karyo(c("46,XX", "46,XX,del(5)(q13q33)", "45,XY,-7"))
risk <- assign_risk(parsed, scheme = "ipssr")
risk[, c("original_karyotype", "ipssr_cyto_risk", "ipssr_cyto_score")]
#> # A tibble: 3 × 3
#>   original_karyotype   ipssr_cyto_risk ipssr_cyto_score
#>   <chr>                <ord>                      <int>
#> 1 46,XX                Good                           1
#> 2 46,XX,del(5)(q13q33) Good                           1
#> 3 45,XY,-7             Poor                           3

# Both schemes, chained
both <- assign_risk(parsed, "ipssr") |> assign_risk("eln2022")
both[, c("original_karyotype", "ipssr_cyto_risk", "eln2022_cyto_risk")]
#> # A tibble: 3 × 3
#>   original_karyotype   ipssr_cyto_risk eln2022_cyto_risk
#>   <chr>                <ord>           <ord>            
#> 1 46,XX                Good            Intermediate     
#> 2 46,XX,del(5)(q13q33) Good            Adverse          
#> 3 45,XY,-7             Poor            Adverse          
```
