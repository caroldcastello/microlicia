# Ideias dos documentos: o que cabe fazer

## Já implementado no Piloto 1.1

- Filtrar ocorrências e estames para as espécies de *Microlicia* com pontas na MCC, guardar a correspondência de nomes e documentar exclusões.
- Checar a árvore MCC, sua ultrametricidade e altura numérica.
- Gerar uma figura da MCC completa, com os 153 táxons nomeados; *Microlicia* fica destacada.
- Resumir os três caracteres categóricos de estames e reconstruir estados ancestrais como uma exploração inicial.
- Cruzar caracteres de estames e resumos geográficos descritivos.
- Gerar hashes, versões do software, tabelas, gráficos e relatório.
- Criar uma etapa opcional de sinalização de coordenadas repetidas por espécie com RuHere, sem exclusão automática.

## Candidatas após revisão dos dados

| Ideia dos documentos | Situação | O que precisa antes |
|---|---|---|
| Comparar a idade da linhagem com estudos prévios | Possível, parcialmente contextualizada | Conferir o XML e escala Ma com o arquivo/execução exata do BEAST; depois comparar os mesmos nós e intervalos de incerteza. |
| BioGeoBEARS: distribuição ancestral e dispersão/vicariância | Ainda não estimada | Definir áreas biogeográficas, completar/validar presença por espécie, revisar pontos e incluir áreas para táxons sem registros; justificar o modelo e limites de áreas. |
| BAMM: shifts de diversificação em Microlicieae ou Microlicia | Não estimada | Confirmar se a MCC é de Microlicia ou da tribo completa, cobertura taxonômica, fração de amostragem e priors/configuração. A MCC sozinha não equivale a uma análise posterior de shifts. |
| MuSSE: diferenças de diversificação entre 3–4 áreas | Não estimada | Codificar estados geográficos definidos, obter frações de amostragem por área e avaliar se o número de espécies/estados suporta os parâmetros. O viés de coleta e a ausência de alguns registros não se corrigem apenas pelo nome do método. |
| Separar regiões do campo rupestre em “cradles”, “museums”, “sources” e “sinks” | Não estimada | Primeiro definir operacionalmente as áreas e a inclusão de espécies de ampla distribuição; a interpretação depende de resultados de diversificação e biogeografia. |
| Diversidade de Microlicieae e de Microlicia para amostragem de Versiane (2019) | Não estimada | Confirmar lista taxonômica, denominadores de riqueza aceita, correspondência entre nomes e cobertura da MCC. |
| Morfometria/ecogeografia do moRphoTaxa | Não executada | O pipeline é voltado a matrizes morfométricas e análises ecogeográficas. O arquivo atual contém somente três caracteres categóricos resumidos por espécie, não medidas contínuas ou dados por espécime; é necessário obter/definir uma matriz morfométrica adequada antes de aplicar PCA/NMDS/modelos de nicho. |

## Uso dos recursos citados

O pipeline moRphoTaxa inspira a organização reprodutível em etapas: harmonização, matriz de caracteres, análises exploratórias e exportação de metadados. Nem todo módulo é transferível ao arquivo atual. As etapas de dados de herbário e limpeza ampla não foram repetidas porque já há uma tabela padronizada; o piloto usa somente a checagem e o filtro de nomes necessários para a MCC.

RuHere foi incluído como etapa opcional de qualidade espacial, começando por duplicatas de coordenadas dentro de espécie. Como a tabela disponível não tem município, estado, país ou classificação por especialista, não executamos validações de país/estado, Flora e Funga do Brasil ou distribuição especializada. As marcações de duplicata são triagem, não evidência de erro.

## Próximas decisões científicas

1. Conferir com o XML efetivamente executado por que a coluna/atributo de unidades está como `years` quando os valores vêm de uma tabela em Ma.
2. Revisar os registros sem coordenadas válidas e os sinalizados por RuHere na versão atualizada das ocorrências; a contagem deve ser regenerada após a troca do arquivo.
3. Decidir áreas e codificação de distribuição antes de BioGeoBEARS, BAMM ou MuSSE.
4. Comparar reconstruções de estames em modelos alternativos e na amostra posterior de árvores antes de usar as probabilidades como resultado principal.
