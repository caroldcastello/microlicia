# Input files

Save the files used by the first analyses here:

- `teste_datacao_MCC.tree` — TreeAnnotator maximum-clade-credibility tree;
- `ocorrencias.csv` — occurrence records;
- `estames.csv` — stamen morphology measurements or coded traits.

Update file paths and exact column names in `config/analysis.yml`. The validation script compares taxon labels literally and writes mismatch reports; it does not rename or edit data.

This repository is public. Data committed here will be publicly accessible. Check data-use terms and coordinate sensitivity before uploading occurrence records. Keep the ~700 MB resampled posterior sample and original 7.64 GB tree logs outside Git; record their generation settings in `data/external/README.md`.
