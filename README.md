# Dashboard Geoquímico — Taió

**Versão 1.0** · Projeto *Taió Plumbing System* (Mestrado em Geologia, UFSC)

🔗 **Acesso:** <https://thubaraleii.github.io/TAIO-DASH/dashboard_geoquimico.html>

Painel interativo de página única (Plotly + Leaflet) com os dados geoquímicos e de campo do sistema intrusivo de Taió,
organizado em três colunas: lista de dados, gráficos e mapa.

## Funcionalidades

**Coluna esquerda**
- Lista unificada de pontos de campo e amostras de geoquímica, com busca por nome ou litologia e seleção que sincroniza gráficos e mapa.
- Resumo dos corpos intrusivos (soleira e dique): área, volume estimado, número de polígonos e classificação Alto/Baixo-Ti.
- Apanhado estatístico de óxidos e elementos-traço (mínimo, máximo, média, mediana, médias por grupo Ti).

**Coluna central**
- **Coluna estratigráfica estimada (furo virtual)** para o ponto selecionado, a partir do modelo estratigráfico da seção 2D.
- **Diagramas gerais:** TAS, AFM e índice de Shand (A/NK × A/CNK).
- **Critérios de Peate et al. (1992):** 12 gráficos (SiO₂, TiO₂, P₂O₅, Fe₂O₃(t), Sr, Ba, Zr, Ti/Zr, Ti/Y, Zr/Y, Sr/Y, Ba/Y), dois por linha, com as janelas dos seis magma-tipos (Urubici, Pitanga, Paranapanema, Ribeira, Esmeralda, Gramado) sobrepostas às amostras.
- **Matrizes de conformidade** com os seis critérios, por grupo (Alto-Ti e Baixo-Ti): pontuação 0–12 por amostra e magma-tipo.
- **Padrões de ETR** normalizados ao manto primitivo (McDonough & Sun, 1995), apenas para amostras com os 14 elementos.
- **Discriminação tectônica Ti/V** (Shervais, 1982) e **diagramas de Harker** (SiO₂ × óxidos maiores).
- Distribuição por litologia e classificação de Ti.

**Coluna direita**
- Mapa Leaflet (Escuro, Rico, Satélite, Relevo, OSM) com camadas de campo, geoquímica e mapa geológico; clique nos pontos sincroniza com a lista e os gráficos.

**Geral:** tema claro/escuro; cor laranja para Alto-Ti e azul para Baixo-Ti em todos os painéis.

## Dados e métodos

- Geoquímica: 41 amostras (rocha total) de afloramentos e furos, com óxidos maiores e elementos-traço; campos derivados calculados na geração:
  Fe₂O₃(t) = FeO(t) × 1,1113; Ti = TiO₂ × 5990; razões Ti/Zr, Ti/Y, Zr/Y, Sr/Y e Ba/Y.
- Critérios e janelas de Peate et al. (1992) aplicados sem alteração; as amostras não são reclassificadas — o grupo Alto/Baixo-Ti vem da planilha.
- Pontuação de conformidade: 1 ponto por critério (de 12) dentro da janela do magma-tipo.

## Uso

Abra `dashboard_geoquimico.html` no navegador (o mapa usa tiles online) ou acesse a versão publicada.

## Regenerar os produtos

O repositório contém o HTML publicado e o script gerador. Para regenerar é necessário o conjunto de dados do projeto
(modelo digital de elevação, litologia, pontos de campo, lineamentos, poços), que **não é distribuído** neste repositório.

```bash
python gerar_dashboard_geoquimico.py
```

Dependências Python: `numpy`, `pandas`, `geopandas`, `shapely`, `pyproj`, `rasterio`, `scipy`, `plotly`, `openpyxl`.

## Estrutura

```
dashboard_geoquimico.html        produto gerado
gerar_dashboard_geoquimico.py    gerador
assets/                          logotipo e ícones
```

## Referências

1. Peate, D. W.; Hawkesworth, C. J.; Mantovani, M. S. M. (1992). Chemical stratigraphy of the Paraná lavas (South America): classification of magma types and their spatial distribution. *Bulletin of Volcanology*, 55, 119–139.
2. Bellieni, G.; Comin-Chiaramonti, P.; Marques, L. S.; Melfi, A. J.; Piccirillo, E. M.; Nardy, A. J. R.; Roisenberg, A. (1984). High- and low-TiO₂ flood basalts from the Paraná plateau (Brazil): petrology and geochemical aspects bearing on their mantle origin. *Neues Jahrbuch für Mineralogie, Abhandlungen*, 150, 273–306.
3. Le Bas, M. J.; Le Maitre, R. W.; Streckeisen, A.; Zanettin, B. (1986). A chemical classification of volcanic rocks based on the total alkali–silica diagram. *Journal of Petrology*, 27, 745–750.
4. Irvine, T. N.; Baragar, W. R. A. (1971). A guide to the chemical classification of the common volcanic rocks. *Canadian Journal of Earth Sciences*, 8, 523–548.
5. Shand, S. J. (1943). *Eruptive Rocks: Their Genesis, Composition, Classification, and Their Relation to Ore-Deposits*. 2nd ed. New York: John Wiley & Sons.
6. Maniar, P. D.; Piccoli, P. M. (1989). Tectonic discrimination of granitoids. *Geological Society of America Bulletin*, 101, 635–643.
7. McDonough, W. F.; Sun, S.-s. (1995). The composition of the Earth. *Chemical Geology*, 120, 223–253.
8. Shervais, J. W. (1982). Ti–V plots and the petrogenesis of modern and ophiolitic lavas. *Earth and Planetary Science Letters*, 59, 101–118.
9. Harker, A. (1909). *The Natural History of Igneous Rocks*. London: Methuen.
10. Milani, E. J.; Melo, J. H. G.; Souza, P. A.; Fernandes, L. A.; França, A. B. (2007). Bacia do Paraná. *Boletim de Geociências da Petrobras*, 15(2), 265–287.
11. Schneider, R. L.; Mühlmann, H.; Tommasi, E.; Medeiros, R. A.; Daemon, R. F.; Nogueira, A. A. (1974). Revisão estratigráfica da Bacia do Paraná. In: *Anais do 28º Congresso Brasileiro de Geologia*, Porto Alegre, SBG, 1, 41–65.
12. Frank, H. T.; Gomes, M. E. B.; Formoso, M. L. L. (2009). Review of the areal extent and the volume of the Serra Geral Formation, Paraná Basin, South America. *Pesquisas em Geociências*, 36(1), 49–57.
13. Peate, D. W. (1997). The Paraná–Etendeka Province. In: Mahoney, J. J.; Coffin, M. F. (eds.), *Large Igneous Provinces*. AGU Geophysical Monograph 100, 217–245.
14. Agafonkin, V. **Leaflet** — an open-source JavaScript library for interactive maps. <https://leafletjs.com>
15. Plotly Technologies Inc. (2015). **Collaborative data science** (Plotly.js). Montréal. <https://plot.ly>
16. Mapzen / Amazon Web Services. **Terrain Tiles** (formato Terrarium). <https://registry.opendata.aws/terrain-tiles/>
17. OpenStreetMap contributors. **OpenStreetMap**. <https://www.openstreetmap.org/copyright>
18. CPRM/SGB — Serviço Geológico do Brasil. **Geologia e Recursos Minerais do Estado de Santa Catarina — Mapa Geológico 1:500.000 (SIG SC)**. <https://geosgb.sgb.gov.br>
19. ANP — Agência Nacional do Petróleo, Gás Natural e Biocombustíveis. **Dados abertos — poços**. <https://www.gov.br/anp>
20. Harris, C. R. et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362.
21. Virtanen, P. et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.
22. Jordahl, K. et al. **GeoPandas**. Zenodo. <https://geopandas.org>; Gillies, S. et al. **Shapely** e **Rasterio**. <https://github.com/shapely/shapely>, <https://github.com/rasterio/rasterio>.

## Licença e uso

- **Código** (scripts de geração e JavaScript das páginas): licença **MIT** — ver [`LICENSE`](LICENSE).
- **Produtos gerados** (HTML, figuras, tabelas derivadas): **CC BY 4.0** — uso, adaptação e redistribuição permitidos com atribuição ao autor.
- **Dados de terceiros** incorporados mantêm as licenças e os termos das fontes originais (tabela em *Créditos e fontes de dados*).
- Os modelos são **interpretativos** e dependem da qualidade e da densidade dos dados de entrada; **não substituem** estudo geológico, hidrogeológico ou geotécnico de detalhe, nem devem ser usados como base única para decisões de engenharia, perfuração ou licenciamento.
- Fornecido “como está”, sem garantias de qualquer tipo.

## Créditos e fontes de dados

| Fonte | Uso | Termos |
|---|---|---|
| Levantamento de campo do autor (pontos, estruturas, geoquímica, MDE de drone) | Dados primários da área de Taió | CC BY 4.0 (autor) |
| CPRM/SGB — Mapa Geológico do Estado de Santa Catarina 1:500.000 (SIG SC) | Geologia regional, falhas, fraturas e diques regionais | Dados abertos CPRM/SGB — citar a fonte |
| ANP — Tabela de poços (dados abertos) | Poços de petróleo/gás (posição e profundidade total) | Dados abertos ANP — citar a fonte |
| Cadastro de poços do Sistema Aquífero Guarani (SAG) | Poços do SAG (posição e profundidade), arquivo KML do projeto | Citar a fonte original do cadastro |
| OpenStreetMap contributors | Rios, estradas e localidades | ODbL 1.0 — © OpenStreetMap contributors |
| Esri World Imagery | Imagem de satélite (tiles e textura) | Termos da Esri — uso com atribuição |
| CARTO / OpenTopoMap | Mapas base do webmap | Atribuição nos próprios mapas (OSM, SRTM, OpenTopoMap CC-BY-SA) |
| AWS Terrain Tiles (formato Terrarium) | Relevo e batimetria regionais (derivados de SRTM, GEBCO e outras fontes) | Atribuição conforme o provedor dos tiles |
| Plotly.js (MIT), Leaflet (BSD-2-Clause) | Bibliotecas de visualização embutidas/carregadas nas páginas | Licenças originais |

## Como citar

> JESUS, A. H. de. **Dashboard Geoquímico — Taió**, versão 1.0. 2026. Projeto *Taió Plumbing System* — Mestrado em Geologia, UFSC. Disponível em: <https://thubaraleii.github.io/TAIO-DASH/dashboard_geoquimico.html>.

## Autoria

**Afonso Henrique de Jesus** — Mestrado em Geologia, Universidade Federal de Santa Catarina (UFSC).
Projeto: *Architecture and Geochemical Evolution of the Taió Intrusive Plumbing System* (Taió, SC — Bacia do Paraná).
