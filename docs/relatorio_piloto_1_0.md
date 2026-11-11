# Microlicia — Piloto de análises 1.0

## Dados e execução

Este relatório é baseado no capítulo “Historical Biogeography of the Melastomataceae”, de Reginato et al. (2022), especialmente Tabela 3 (p. 94), nas ocorrências e estames fornecidos e na árvore MCC. O workflow mantém dados de ocorrência somente para táxons encontrados na MCC.

Análises feitas com a árvore MCC, ocorrências previamente filtradas para táxons da árvore e matriz categórica de estames. O script `scripts/pilot_1_0.py` gera tabelas, figuras e metadados. O XML do BEAST não é versionado; seus rótulos de unidade e priors foram registrados em `data/metadata/fonte_datação.json`.

## Resultados de validação

- Árvore: 153 pontas únicas; ultramétrica (amplitude raiz–ponta = 1.42e-14); altura da raiz = 28.048261 unidades numéricas da árvore. O XML local declara units=['years']; médias/priores avaliados: 29.47, 30.03, 22.39.
- [Inferência] A escala pretendida é Ma: 29,47 (Marcetieae), 30,03 (Rhexieae) e 22,39 (Lavoisiereae) correspondem às idades da Tabela 3 do capítulo de Reginato et al. (2022), cuja coluna está em Ma. A altura da raiz da MCC fica, portanto, provisoriamente interpretada como 28,048 Ma. O XML, porém, declara `years`; essa divergência deve ser reconciliada com a configuração efetivamente executada no BEAST antes de publicação ou reexecução da datação.
- Estames: 153 táxons nomeados; 153 correspondem exatamente às pontas da árvore. Os três caracteres têm dados para todos os táxons nomeados.
- Ocorrências: 10208 registros; 10199 pares de coordenadas completos dentro dos limites numéricos de latitude/longitude; 9 sem par válido.
- Nomes: a análise recebeu ocorrências já filtradas para as pontas da MCC após remover apenas espaços entre `m_` e o epíteto. Foram mantidos 10208 registros de 142 espécies. Antes da filtragem, 24 táxons (117 registros) não tiveram correspondência; a lista está em `data/metadata/taxa_ocorrencias_excluidos_sem_correspondencia.csv` e nos resultados `occurrence_taxa_unmatched.csv`.
- Cobertura espacial ligada à árvore: 10208 registros de 142 táxons; 10199 registros têm coordenadas válidas.

## Estames e reconstrução ancestral

As frequências observadas estão em `stamen_trait_frequencies_and_root_posterior.csv`. A reconstrução ancestral foi feita por máxima verossimilhança sob um modelo de transições iguais (ER), na MCC, com probabilidades iguais para os estados na raiz antes de observar os dados. As probabilidades por nó estão nos três arquivos `ancestral_reconstruction_trait*.csv`.

- Probabilidade posterior estimada na raiz para **dimorfismo estaminal**: dimórficos = 0.955; subisomórficos = 0.013; isomórficos = 0.032.
- Probabilidade posterior estimada na raiz para **padrão de coloração dos estames**: concolores = 0.544; bicolores = 0.456.
- Probabilidade posterior estimada na raiz para **esporângios das anteras**: tetrasporangiadas = 0.918; poliesporangiadas = 0.082.
[Inferência] Como a escala pretendida é Ma, as taxas `q` podem ser lidas por Ma apenas se a divergência entre a fonte de calibração e o rótulo `years` do XML for resolvida. A reconstrução depende do modelo e não incorpora variação entre árvores.

## Ocorrências e associação descritiva

`occurrence_species_summary_tree_matched.csv` resume número de registros, coordenadas únicas, centroides aritméticos e amplitude em graus por táxon coincidente com a árvore. `stamen_traits_with_occurrence_summaries.csv` cruza os caracteres e esses resumos. `exploratory_trait_geography_descriptives.csv` apresenta medianas e amplitudes de latitude/longitude e esforço de coleta por estado morfológico; são estatísticas descritivas, sem teste causal. O critério de coordenada válida nesta etapa é apenas numérico (latitude entre −90 e 90, longitude entre −180 e 180); isso não confirma a localização geográfica de cada registro.

## Figuras

- `occurrence_coordinates_scatter.png`: pontos em longitude/latitude decimal, sem projeção ou limites políticos.
- `stamen_trait_state_counts.png`: frequências dos estados dos três caracteres.
- `exploratory_trait_latitude.png`: distribuição descritiva da latitude média por estado nos dois primeiros caracteres.

## Limites do Piloto 1.0

- Os 24 nomes sem correspondência (117 registros) foram excluídos das ocorrências de entrada após normalização apenas de espaços em `m_`; a tabela de exclusão preserva nomes e contagens, sem inferir sinônimos ou corrigir grafias.
- Nove registros das espécies mantidas não têm par de coordenadas válido pelos limites numéricos e foram excluídos dos resumos espaciais; precisam de revisão.
- Não há CRS documentado, áreas biogeográficas codificadas nem camadas ambientais neste piloto. Não foram estimadas áreas em km², modelos ambientais ou biogeografia histórica. A execução usa a MCC; a amostra posterior não foi processada nesta versão.
- A reconstrução ER é a análise exploratória inicial; modelos alternativos e sensibilidade precisam ser avaliados antes de inferência final, principalmente para o caráter com um único táxon no estado “isomórficos”.

## Ideias dos documentos ainda não executadas

BioGeoBEARS, BAMM e MuSSE não foram executados neste piloto. Faltam a codificação de áreas por espécie, o tratamento de cobertura amostral e a confirmação do escopo da MCC para sustentar essas análises. A proposta BAMM de um shift compartilhado pelo campo rupestre é uma hipótese a testar, não uma conclusão deste piloto. O moRphoTaxa pipeline orienta a organização de dados e análises morfométricas/ecogeográficas, mas as três variáveis atuais são caracteres categóricos já resumidos por espécie; não há matriz de medidas morfométricas por espécime para rodar PCA/NMDS ou modelos de nicho. RuHere foi incluído como triagem opcional de coordenadas duplicadas; as marcações não removem linhas.

## Arquivos de saída

`data/metadata/pilot_1_0_metadata.json` registra hash SHA-256 e tamanho dos arquivos usados, parâmetros gerais e limitações. A execução no GitHub Actions gera e disponibiliza as tabelas, figuras, relatório e manifesto como artefato. A tabela de entrada contém apenas ocorrências das espécies correspondentes à MCC.
