# Workflow: dated tree, stamen morphology, and occurrences

## Scope

Start from the existing dated MCC tree; do not rerun BEAST for this initial workflow. The MCC tree supports single-tree analyses. Keep the resampled posterior tree file for later analyses that should propagate tree uncertainty.

## Stage 1 — organize and validate inputs

| Input | Minimum content |
|---|---|
| `data/input/teste_datacao_MCC.tree` | TreeAnnotator MCC tree with branch lengths |
| `data/input/ocorrencias.csv` | One row per occurrence, taxon name, longitude, latitude |
| `data/input/estames.csv` | Taxon name and measured or coded stamen traits |

Set actual file and column names in `config/analysis.yml`. Taxon matching is literal: the script reports mismatches and does not silently standardize names.

Before statistical analysis, document taxonomic spelling/synonym decisions and sources; units and definitions for each quantitative trait; definitions of categorical states; whether rows represent species summaries, specimens, or repeated measurements; and sample size per species where available. Do not automatically average repeated measurements: choose and document a species-level summary or a hierarchical approach after reviewing the sampling design.

## Stage 2 — tree checks

Check that the MCC tree can be read, has unique tip labels and branch lengths, and is ultrametric within numerical tolerance. Inspect root age and units against the BEAST configuration before reporting dates in years or millions of years. A successful file read does not validate the calibration or dating model.

For the posterior sample, retain a metadata note stating the original logging interval, burn-in removed, resampling interval, and resulting tree count. The filename alone does not establish these settings.

## Stage 3 — stamen traits on the MCC tree

For quantitative traits: inspect ranges, missingness, and distributions; calculate a species-level value appropriate to the sampling design; assess phylogenetic signal (for example, Pagel's lambda); visualize trait distribution on the tree and in trait space; and reconstruct ancestral values only when coverage and assumptions are adequate.

For categorical traits: define biologically meaningful states; report state frequencies and missing taxa; compare plausible transition models (for example, equal-rates and all-rates-different); and report ancestral-state uncertainty rather than a single state as certain.

The actual models and code will be fixed after inspecting the trait table and sampling structure. A reconstructed association alone does not demonstrate adaptation.

## Stage 4 — occurrence data

First summarize records per species, coordinate completeness, duplicate coordinates, and geographic coverage. Map records only after confirming the coordinate reference system and cleaning rules.

Potential next analyses depend on the fields and question: species range size/geographic spread; comparison of stamen traits across pre-defined regions or environments; or historical biogeography using a justified area coding scheme and adequate taxon sampling. Select a model only after the geographic areas, outgroup/taxon sampling, and occurrence-cleaning choices are specified.

## Stage 5 — combine traits, geography, and phylogeny

After Stages 1–4, formulate specific hypotheses and select tests from the actual columns and sampling design. A phylogenetic comparative model may test association between a quantitative stamen trait and a defined geographic variable. Correlation does not by itself establish causal adaptation.

## Stage 6 — propagate tree uncertainty

Repeat selected estimates across the resampled posterior trees and summarize variation among trees. Label MCC results as single-tree analyses. Keep the approximately 700 MB posterior sample outside Git history and use an external storage or compute location when needed.

## Outputs and reproducibility

Keep raw inputs unchanged. Store summaries, figures, and model outputs under `results/`. Record software versions, configuration, exclusions, and taxonomic decisions. Do not publish coordinates with access restrictions or sensitive-location concerns.
