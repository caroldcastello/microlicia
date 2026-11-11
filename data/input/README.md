# Entradas do piloto Microlicia 1.1

- `teste_datacao_MCC.tree` — árvore MCC datada com 153 pontas. A figura da árvore usa todos os táxons e identifica os nomes científicos completos; as pontas de *Microlicia* ficam destacadas.
- `ocorrencias.csv` — tabela padronizada de ocorrências. As análises espaciais usam apenas registros cujo táxon corresponde a uma ponta de *Microlicia*. Registros de outros gêneros e nomes sem correspondência são listados em arquivos de auditoria e não entram nas análises focais.
- `estames.csv` — matriz categórica de estames. As frequências e reconstruções usam apenas táxons de *Microlicia* presentes na MCC.

A MCC original não é alterada. O script corrige em memória os rótulos `m_bellineloi` → `m_belinelloi` e `m_diamantinensis` → `m_diamantinesis`, conforme os nomes corrigidos nas ocorrências.

A amostra posterior (~700 MB) e os logs do BEAST não são necessários para gerar a figura da MCC. Coordenadas de ocorrência estão versionadas neste repositório público; mantenha os dados aqui somente se isso estiver de acordo com os termos de uso e eventuais restrições de localização.
