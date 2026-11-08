# Microlicia: dated-tree analyses

Reproducible workflow for analyses combining a dated phylogeny, species occurrence records, and stamen morphology.

## Where to save the inputs

- `data/input/`: MCC tree, occurrence table, and stamen-trait table used by the first analyses.
- `data/external/`: notes and pointers for large files kept outside Git, including the resampled posterior tree sample and original BEAST tree logs.
- `results/`: generated summaries, figures, and model outputs.

**This repository is public. Any data committed here will also be public.** Check data-use terms and coordinate sensitivity before committing occurrence records. The resampled tree sample (~700 MB) and original 7.64 GB logs should remain outside the repository.

## Start here

1. Read [the workflow](docs/workflow.md).
2. Place the MCC tree and the two CSV files in `data/input/`.
3. Set actual file and column names in `config/analysis.yml`.
4. Install R packages `ape` and `yaml`, then run from the repository root:
   ```r
   install.packages(c("ape", "yaml"))
   source("scripts/00_validate_inputs.R")
   ```
5. Review the reports written to `results/00_input_checks/` and resolve taxon-name mismatches before statistical analyses.

The checks do not rename taxa, infer missing values, or modify source data. Statistical models will be configured after inspecting the actual traits, replication structure, and occurrence fields.
