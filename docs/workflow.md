# Workflow do piloto 1.0

## 1. Filtro de espécies

As ocorrências foram cruzadas com as 153 pontas da MCC. A regra de correspondência remove somente espaços entre `m_` e o epíteto. Foram mantidos 10.208 registros de 142 espécies. Os 117 registros de 24 nomes sem correspondência foram excluídos da tabela de análise, sem sinonimização ou correção de grafia. A tabela de auditoria de nomes preserva as contagens excluídas.

A tabela de estames tem um registro para cada uma das 153 pontas e três caracteres categóricos completos. A coluna de táxon é o primeiro campo do CSV e veio sem nome de cabeçalho.

## 2. Verificações da MCC e datas

O script confirma ponta única, comprimentos de ramo e ultrametricidade. A altura da raiz foi calculada em unidades numéricas da árvore. A Tabela 3 do capítulo de Reginato et al. (2022) identifica em Ma as idades 29,47 Ma (coroa de Marcetieae), 30,03 Ma (Rhexieae) e 22,39 Ma (Lavoisiereae), coincidentes com os valores dos priors no XML disponível.

[Inferência] Isso aponta que a escala pretendida é Ma e que a raiz da MCC representa aproximadamente 28,048 Ma. O XML traz `units="years"`; o valor da tabela e o rótulo do XML divergem. O relatório conserva essa ressalva e não altera o XML nem reexecuta o BEAST.

## 3. Estames

O script calcula frequências e reconstrói estados ancestrais por máxima verossimilhança sob CTMC de taxas iguais (ER), com priori uniforme na raiz, usando apenas a MCC. As probabilidades por nó são estimativas dependentes do modelo. Antes de resultados finais, comparar ER com modelos alternativos e resumir a incerteza entre árvores posteriores.

## 4. Ocorrências

Os resumos usam registros com latitude e longitude numéricas dentro dos limites globais. O filtro não valida que uma coordenada esteja no lugar certo. Os nove pares fora dos limites ou incompletos são listados para revisão. Resumos em graus são descritivos; não estimam área geográfica e não substituem revisão espacial.

O script opcional `scripts/ruhere_duplicate_qc.R` usa `RuHere::flag_duplicates()` para sinalizar coordenadas repetidas dentro de espécie. Ele exporta os sinalizadores e mantém todas as linhas: a decisão de remover qualquer registro exige inspeção de vouchers e metadados.

## 5. Reprodutibilidade

`python scripts/pilot_1_0.py` gera tabelas e figuras em `results/`; a execução no GitHub Actions salva a pasta como artefato. O manifesto guarda versões, parâmetros, hashes dos arquivos e exclusões. A árvore posterior e os logs extensos não são necessários para o piloto MCC e permanecem fora do Git.
