# Taió Plumbing System — Dashboard Geoquímico

Painel geoquímico (Plotly + Leaflet, HTML autônomo, roda offline sem
servidor) do sistema intrusivo de Taió — sill/soleira e dique de diabásio
intrudindo formações sedimentares da Bacia do Paraná.

Lista/busca dos pontos de campo sincronizada com um mapa Leaflet (basemaps
satélite/relevo-topográfico/mapa geológico real CPRM) e 11 diagramas
geoquímicos de classificação (TAS, AFM, Shand's Index, MgO×TiO₂,
Ti/Y×Ti/Zr, Fe₂O₃×TiO₂, Zr/Y×Sr, Zr/Y×Ti/Zr, Sr×Ti/Y, litologia e Ti Taió),
com campos de referência digitalizados a partir da literatura (Fontoura,
TCC 2024; Peate et al. 1992/1997) — Taió ainda não tem óxido/traço bruto
próprio pra plotar, então os diagramas mostram só os campos de referência.

Abra `dashboard_geoquimico.html` direto no navegador.

Complemento do modelo 3D: [taio-plumbing-system-3d](https://github.com/Thubaraleii/taio-plumbing-system-3d).

Gerado por `gerar_dashboard_geoquimico.py` (Python, Plotly/GeoPandas) — o
script depende dos dados do projeto completo (catálogo de pontos de campo,
mapa geológico CPRM) e não roda de forma standalone fora dessa estrutura;
está incluído aqui só como referência/histórico do código.

Criado por Afonso Henrique de Jesus.
