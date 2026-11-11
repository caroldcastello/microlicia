# Workflow Microlicia — Piloto 1.1

## Escopo

As análises de ocorrências, morfologia de estames e reconstrução de estados ancestrais são restritas às pontas de *Microlicia* da MCC. A figura da árvore datada mostra todas as 153 pontas da árvore, incluindo os gêneros usados como grupos externos na datação. Os dois rótulos corrigidos nas ocorrências são recodificados em memória na árvore e na matriz de estames, sem modificar o arquivo MCC original.

## Entradas

- `data/input/teste_datacao_MCC.tree`
- `data/input/ocorrencias.xlsx.csv`
- `data/input/estames.csv`
- `data/metadata/fonte_datação.json` (referência às unidades/priores quando o XML do BEAST não está no repositório)

## Execução

```bash
python scripts/pilot_1_0.py
```

Para gerar somente a figura da MCC completa:

```bash
python scripts/pilot_1_0.py --tree-figure-only
```

O workflow do GitHub Actions salva a pasta `results/` como artefato. As análises comparativas usam a subárvore induzida pelas pontas de *Microlicia*; a figura da MCC usa a árvore original e todos os rótulos.

## Figuras

- `dated_tree_all_taxa.png` e `dated_tree_all_taxa.pdf`: árvore datada com 153 nomes científicos completos; *Microlicia* destacada em verde. O eixo usa unidades numéricas da árvore até que a escala pretendida seja reconciliada com a configuração BEAST.
- `occurrence_coordinates_scatter.png`: registros de *Microlicia* com coordenadas numericamente válidas.
- `stamen_trait_state_counts.png`: frequências dos três caracteres de estames em *Microlicia*.
- `exploratory_trait_latitude.png`: resumo exploratório de latitude por estados morfológicos de *Microlicia*.

## Auditoria

`occurrence_taxon_crosswalk.csv` separa correspondências com a MCC completa e inclusão nas análises de *Microlicia*. `occurrence_outgroups_excluded_from_analyses.csv` lista registros dos demais gêneros que foram excluídos do foco. `occurrence_taxa_unmatched.csv` lista nomes não encontrados na árvore, e `microlicia_taxa_without_occurrences.csv` lista pontas de *Microlicia* sem registro correspondente.
