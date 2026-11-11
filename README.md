# Microlicia — piloto 1.0

Repositório de análises reproduzíveis com a árvore datada MCC, caracteres de estames e dados de ocorrência filtrados para os táxons amostrados na árvore.

## Conteúdo de dados

`data/input/` contém a MCC, a tabela de estames e a tabela de ocorrências após filtro taxonômico. A coluna `taxon_tree` guarda o nome comparável à ponta da árvore. A única harmonização automática aplicada aos nomes foi remover espaços entre `m_` e o epíteto; não foram atribuídos sinônimos nem corrigidas grafias.

O arquivo de ocorrências contém 10.208 registros de 142 táxons da MCC. A MCC tem 153 pontas; nove espécies não têm registro de ocorrência no conjunto filtrado. Foram deixados de fora 117 registros de 24 nomes sem correspondência; a auditoria dos nomes e contagens está em `data/metadata/taxa_ocorrencias_excluidos_sem_correspondencia.csv`.

**Este repositório é público.** Os arquivos em `data/input/` ficam publicamente acessíveis. A tabela contém coordenadas de ocorrências; confira termos de uso e eventuais restrições de localização antes de manter os dados aqui.

## Executar o piloto

O script Python principal gera tabelas, gráficos, relatório e metadados em `results/`:

```bash
python -m pip install -r requirements.txt
python scripts/pilot_1_0.py
```

Para refazer o processamento em outra pasta, defina `PILOT_INPUT_DIR` e `PILOT_OUTPUT_DIR`. A execução requer os arquivos em `data/input/` e o registro de calibração em `data/metadata/`.

Também é possível executar o workflow no GitHub em **Actions → Microlicia pilot 1.0 → Run workflow**. O relatório e os arquivos de resultados serão disponibilizados como artefato da execução.

## Resultados e interpretação

- Relatório do piloto: [`docs/relatorio_piloto_1_0.md`](docs/relatorio_piloto_1_0.md).
- Metadados e hashes da execução: [`data/metadata/pilot_1_0_metadata.json`](data/metadata/pilot_1_0_metadata.json).
- Fonte dos valores de calibração e observação sobre escala: [`data/metadata/fonte_datação.json`](data/metadata/fonte_datação.json).
- Próximas análises possíveis e limites: [`docs/proximas_analises.md`](docs/proximas_analises.md).

[Inferência] A escala pretendida das idades é Ma: os três valores de calibração coincidem com as idades da Tabela 3 de Reginato et al. (2022), expressas em Ma. O XML local indica `years`; essa divergência deve ser resolvida com a configuração efetivamente usada no BEAST antes de tratar as idades como confirmadas.

O arquivo de árvores posterior (~700 MB) e os logs grandes do BEAST não estão incluídos. O piloto atual usa somente a MCC e não propaga incerteza entre árvores.
