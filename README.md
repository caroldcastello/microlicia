# Microlicia — piloto de análises 1.3

Repositório de análises reproduzíveis com a árvore datada MCC, dados de ocorrência e caracteres de estames. As análises de ocorrências, estames e reconstrução ancestral têm foco somente em *Microlicia*. A figura da MCC inclui todos os táxons da árvore.

## Conteúdo de dados

`data/input/` contém a MCC, a matriz de estames e o arquivo de ocorrências. O script reconhece `ocorrencias.xlsx.csv` ou `ocorrencias.csv`; se ambas existirem, usa primeiro `ocorrencias.xlsx.csv`. A correspondência pode usar `taxon_tree` quando essa coluna existir. Nomes com espaço entre `m_` e o epíteto são harmonizados, e os dois rótulos indicados no registro do projeto são corrigidos em memória para alinhar árvore, estames e ocorrências.

A MCC tem 153 pontas: 145 de *Microlicia* e oito de outros gêneros. Esses outros gêneros continuam na árvore datada e na figura, mas são excluídos das análises comparativas de *Microlicia*. O script grava auditorias para nomes não encontrados e pontas de *Microlicia* sem ocorrência.

**Este repositório é público.** Os arquivos em `data/input/` ficam publicamente acessíveis. A tabela contém coordenadas de ocorrências; confira termos de uso e eventuais restrições de localização antes de manter os dados aqui.

## Executar o piloto

O script Python principal gera tabelas, figuras, relatório e metadados em `results/`:

```bash
python -m pip install -r requirements.txt
python scripts/pilot_1_0.py
```

Para gerar apenas a árvore datada completa, sem executar as análises de ocorrências e estames:

```bash
python scripts/pilot_1_0.py --tree-figure-only
```

Essa etapa cria `results/dated_tree_all_taxa.png` e `results/dated_tree_all_taxa.pdf`, com nomes científicos completos e *Microlicia* destacada. O eixo temporal mostra as unidades numéricas da MCC até a divergência entre a escala das calibrações e o rótulo do XML ser reconciliada.

Para refazer o processamento em outra pasta, defina `PILOT_INPUT_DIR` e `PILOT_OUTPUT_DIR`. A execução requer os arquivos em `data/input/` e o registro de calibração em `data/metadata/`.

Também é possível executar o workflow no GitHub em **Actions → Microlicia pilot 1.3 → Run workflow**. O relatório e os arquivos de resultados serão disponibilizados como artefato da execução.

## Resultados e interpretação

- O relatório histórico do piloto 1.0 está em [`docs/relatorio_piloto_1_0.md`](docs/relatorio_piloto_1_0.md); o relatório atualizado é gerado em `results/relatorio_piloto_1_3.md` ao executar a versão 1.3.
- Metadados e hashes da execução: [`data/metadata/pilot_1_0_metadata.json`](data/metadata/pilot_1_0_metadata.json).
- Fonte dos valores de calibração e observação sobre escala: [`data/metadata/fonte_datação.json`](data/metadata/fonte_datação.json).
- Próximas análises possíveis e limites: [`docs/proximas_analises.md`](docs/proximas_analises.md).

[Inferência] A escala pretendida das idades é Ma: os três valores de calibração coincidem com as idades da Tabela 3 de Reginato et al. (2022), expressas em Ma. O XML local indica `years`; essa divergência deve ser resolvida com a configuração efetivamente usada no BEAST antes de tratar as idades como confirmadas.

O arquivo de árvores posterior (~700 MB) e os logs grandes do BEAST não estão incluídos. O piloto atual usa somente a MCC e não propaga incerteza entre árvores.
